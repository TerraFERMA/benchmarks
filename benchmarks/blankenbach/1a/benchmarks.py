name = "blankenbach"
ncells = [16, 32, 64, 128, 256]

import os
import subprocess
import sys
import libspud
from buckettools.statfile import parser
from math import sqrt
import hashlib
import shutil
from string import Template as template

tf_repo_path = os.environ["TF_CMAKE_PATH"]

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
    
    tempgeo = file("../../../src/transfinite_square_template.geo", 'r')
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
  table0.append("\\begin{tabular}{c|cccccc}\n")
  table0.append("    & $\Nu$ & $v_{rms}$ & $q_1$ & $q_2$ & $T_{e}$ & $z_{e}$ \\\\\n")
  table0.append("\hline\n")

  libspud.load_options(name+".tfml")
  ioname = libspud.get_option("/io/output_base_name")
  Ra = libspud.get_option("/system::Stokes/coefficient::RayleighNumber/type::Constant/rank::Scalar/value::WholeMesh/constant")
    
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

    nu = -1.0*(stat["Stokes"]["Temperature"]["TopSurfaceIntegral"][-1])
    v_rms = sqrt(stat["Stokes"]["Velocity"]["L2NormSquared"][-1])*Ra
    q1 = -det["TemperatureGradient"]["TemperatureGradient"]["TopLeft"][-1][-1]
    q2 = -det["TemperatureGradient"]["TemperatureGradient"]["TopRight"][-1][-1]
    q3 = -det["TemperatureGradient"]["TemperatureGradient"]["BottomRight"][-1][-1]
    q4 = -det["TemperatureGradient"]["TemperatureGradient"]["BottomLeft"][-1][-1]
    detlen = len(det["Stokes"]["Temperature"]["Array"][:,-1])
    expos1 = det["Stokes"]["Temperature"]["Array"][0:detlen/2,-1].argmin()
    Te1 = det["Stokes"]["Temperature"]["Array"][0:detlen/2,-1][expos1]
    ze1 = det["Array"]["position_1"][0:detlen/2,-1][expos1]
    expos2 = det["Stokes"]["Temperature"]["Array"][detlen/2:,-1].argmax()
    Te2 = det["Stokes"]["Temperature"]["Array"][detlen/2:,-1][expos2]
    ze2 = det["Array"]["position_1"][detlen/2:,-1][expos2]

    table0.append("{0:2d}$\\times${0:2d} & {1:.3f} & {2:.3f} & {3:.3f} & {4:.3f} & {5:.3f} & {6:.3f} \\\\\n".format(n, nu, v_rms, q1, q2, Te1, ze1))

    os.chdir("../")
  
  libspud.clear_options()

  table0.append("\hline\n")
  table0.append("Benchmark & 4.884 & 42.865 & 8.059 & 0.589 & 0.422 & 0.225 \\\\\n")
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

