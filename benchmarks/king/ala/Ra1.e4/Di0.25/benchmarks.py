name = "king"
ncells = [16, 32, 64, 128]

approx = "ALA"
comparisoncodes = ['UM', 'VT', 'CU', 'CT', 'CZ', 'KS']

import os
import subprocess
import sys
import libspud
from buckettools.statfile import parser
from math import sqrt
import hashlib
import shutil
from string import Template as template
import pickle

def generate():
  libspud.load_options(name+".tfml")
  meshbasename = libspud.get_option("/geometry/mesh::Mesh/source::File/file")
  libspud.clear_options()
  
  for n in ncells:
    dir_name = "run_"+`n`
    try:
      os.mkdir(dir_name)
    except OSError:
      pass
    chdir(dir_name)
    
    tempgeo = file("../../../../../src/transfinite_square_template.geo", 'r')
    templines = tempgeo.readlines()
    tempgeo.close()
    
    lines = []
    for templine in templines:
      line = template(templine).substitute({'ncells':n, 'factor':0.01})
      lines.append(line)
    
    geo = file(meshbasename+"_temp.geo", 'w')
    geo.writelines(lines)
    geo.close()
    
    try:
      checksum = hashlib.md5(open(meshbasename+".geo").read()).hexdigest()
    except:
      checksum = None

    if checksum != hashlib.md5(open(meshbasename+"_temp.geo").read()).hexdigest():
      # file has changed
      shutil.copy(meshbasename+"_temp.geo", meshbasename+".geo")
      try:
        subprocess.check_call(["gmsh", "-2", "-algo", "del2d", meshbasename+".geo"])
      except:
        print "ERROR while calling gmsh in directory ", dir_name
        os.chdir("../")
        sys.exit(1)
      
      try:
        subprocess.check_call(["dolfin-convert", meshbasename+".msh", meshbasename+".xml"])
      except:
        print "ERROR while calling dolfin-convert in directory ", dir_name
        os.chdir("../")
        sys.exit(1)
      
      extensions = [".xml", "_facet_region.xml"]
      for ext in extensions:
        try:
          subprocess.check_call(["gzip", "-f", meshbasename+ext])
        except:
          print "ERROR while calling gzip in directory ", dir_name, " on extension ", ext
          os.chdir("../")
          sys.exit(1)
    
    os.chdir("../")

def configure():

  dir_name = "build"
  try:
    subprocess.call(["tfbuild", name+".tfml", "-d", dir_name])
  except:
    print "ERROR while calling tfbuild for  directory ", dir_name
    sys.exit(1)


def build():

  dir_name = "build"
  chdir(dir_name)
  try:
    subprocess.check_call(["make"])
  except:
    print "ERROR while calling make in directory ", dir_name
    os.chdir("../")
    sys.exit(1)
  os.chdir("../")

def run():
  for n in ncells:
    dir_name = "run_"+`n`
    chdir(dir_name)
    try:
      print "running in ",dir_name
      subprocess.check_call(["../build/"+name, "-vINFO", "-l", "../"+name+".tfml"])
    except:
      print "ERROR while calling running in directory ", dir_name
      os.chdir("../")
      sys.exit(1)
    os.chdir("../")
  
def extract():
  table0 = []
  table0.append("\\begin{tabular}{c|ccccccc}\n")
  table0.append("    & $\Nu$ & $v_{rms}$ & $\max(u|_{z=1})$ & $<u|_{z=1}>$ & $<T + \\bar{T}>$ & $<\phi>$ & $<W>$ \\\\\n")
  table0.append("\hline\n")

  libspud.load_options(name+".tfml")
  ioname = libspud.get_option("/io/output_base_name")
  Ra = libspud.get_option("/system::Stokes/coefficient::RayleighNumber/type::Constant/rank::Scalar/value::WholeMesh/constant")
  Di = libspud.get_option("/system::Stokes/coefficient::DissipationNumber/type::Constant/rank::Scalar/value::WholeMesh/constant")
    
  for n in ncells:
    dir_name = "run_"+`n`
    chdir(dir_name)
    
    try:
      stat = parser(ioname+".stat")
    except IOError:
      os.chdir("../")
      continue
    try:
      det  = parser(ioname+".det")
    except IOError:
      os.chdir("../")
      continue

    nu = -1.0*(stat["Stokes"]["Temperature"]["FullTopSurfaceIntegral"][-1])
    v_rms = sqrt(stat["Stokes"]["Velocity"]["L2NormSquared"][-1])*Ra
    v_surf_max = abs(det["Stokes"]["Velocity_0"]["Array"][:,-1]).max()*Ra
    v_surf_int = stat["Stokes"]["Velocity"]["TopSurfaceIntegral"][-1]*Ra
    T_int = stat["Stokes"]["Temperature"]["FullIntegral"][-1]
    phi_int = stat["ViscousDissipation"]["ViscousDissipation"]["Integral"][-1]*Ra*Di
    W_int = stat["WorkDone"]["WorkDone"]["Integral"][-1]*Ra

    table0.append("{0:2d}$\\times${0:2d} & {1:.3f} & {2:.3f} & {3:.3f} & {4:.3f} & {5:.3f} & {6:.3f} & {7:.3f} \\\\\n".format(n, nu, v_rms, v_surf_max, v_surf_int, T_int, phi_int, W_int))

    os.chdir("../")
  
  libspud.clear_options()

  Raind = int(round(Ra))
  Diind = '%d'%Di if Di == int(Di) else '%s'%Di
  f = file("../../../tables.dict", 'r')
  tables = pickle.load(f)
  f.close()
   
  table0.append("\hline\n")
  for code,values in tables[approx][Raind][Diind].iteritems():
    if len(values)>0 and code in comparisoncodes:
      line = "Benchmark ("+code+") & "
      for v in range(len(values)-1):
        line += '%.3f'%values[v]+" & "
      line += '%.3f'%values[-1]+" \\\\\n"
      table0.append(line)
  table0.append("\\end{tabular}\n")

  for line in table0:
    print line[0:-1]

  filename = "table0.tex"
  filehandle = file(filename, 'w')
  filehandle.writelines(table0)
  filehandle.close()

def chdir(dirname):
  try:
    os.chdir(dirname)
  except OSError:
    print "ERROR no such directory "+dirname+"."
    sys.exit(1)

