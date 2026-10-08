"""Whole native BNDINP/IODISK input oracle; STOP before any PDE or AUXVAL return."""
import argparse,hashlib,json,math,os,shutil,subprocess
from pathlib import Path
import numpy as np
from scipy.io import FortranFile

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--forward-input',type=Path);a=p.parse_args();a.source=a.source.resolve();a.output=a.output.resolve();a.output.mkdir(parents=True)
 here=Path(__file__).resolve().parent;build=a.output/'build';build.mkdir()
 for s in a.source.iterdir():
  if s.is_file()and(s.suffix in ['.f','.inc']or s.name=='makefile'):shutil.copy2(s,build/s.name)
 original=(build/'chease.f').read_text();begin=original.index('         SUBROUTINE BNDINP');end=original.index('C*DECK C2SX03',begin)
 block=original[begin:end];assert block.count('         CALL IODISK(33)')==1
 block=block.replace("         INCLUDE 'DECLAR.inc'","         INCLUDE 'DECLAR.inc'\n         INCLUDE 'COMDIM.inc'\n         INCLUDE 'COMBND.inc'\n         INCLUDE 'COMCON.inc'\n         INCLUDE 'COMBAL.inc'\n         INCLUDE 'COMPHY.inc'")
 block=block.replace('         CALL IODISK(33)',"""         CALL IODISK(33)
         NPPF1=NPPF+1
         OPEN(17,FILE='parsed.bin',FORM='UNFORMATTED')
         WRITE(17) NBPS,NWBPS,NDATA,NPPF1,NSTTP
         WRITE(17) ASPCT,RZ0C,PREDGE,
     & RRBPS(1:NBPS,1),RZBPS(1:NBPS,1),
     & FCSM(1:NPPF1),RPPF(1:NPPF1),RFUN(1:NPPF1)
         CLOSE(17)
         STOP 'NATIVE INPUT ORACLE: NO PDE'
""")
 (build/'chease.f').write_text(original[:begin]+block+original[end:])
 flags='-O0 -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-none -std=legacy -fallow-argument-mismatch -mcmodel=medium'
 command=['make','-j1','F95=gfortran','F95FLAGS='+flags]
 with(build/'build.log').open('wb')as f:r=subprocess.run(command,cwd=build,stdout=f,stderr=subprocess.STDOUT)
 assert r.returncode==0,'Native Make failed: '+str(build/'build.log')
 subprocess.run(['python',str(here/'generate_inputs.py'),str(a.output/'generated'),'--variant','mars'],check=True,capture_output=True)
 cases={}
 for shape,kappa in [('circular',1.),('shaped',1.7)]:
  for selector in [4,1]:
   d=a.output/(shape+'_selector'+str(selector));d.mkdir();template=a.output/'generated'/shape
   lines=(template/'EXPEQ').read_text().splitlines();assert lines[3]=='513 1 1' and lines[517]=='101 4'
   lines[517]='101 '+str(selector)
   if selector==1:lines[720:821]=['-1.2500000000000000e-01']*101
   (d/'EXPEQ').write_text('\n'.join(lines)+'\n');(d/'chease_namelist').write_text((template/'chease_namelist').read_text().replace('NSTTP=4','NSTTP='+str(selector)))
   cases[d.name]=(d,kappa,selector,None)
 if a.forward_input:
  d=a.output/'actual_forward';d.mkdir()
  for name in ['EXPEQ','chease_namelist']:shutil.copyfile(a.forward_input/name,d/name)
  cases[d.name]=(d,1.,1,np.array([float(s)for s in(d/'EXPEQ').read_text().splitlines()[720:821]]))
 for defect in ['missing_wall_header','missing_title_records']:
  d=a.output/defect;d.mkdir();template=a.output/'generated'/'circular'
  lines=(template/'EXPEQ').read_text().splitlines()
  deck=(template/'chease_namelist').read_text()
  if defect=='missing_wall_header':lines[3]='513'
  else:deck='\n'.join(deck.splitlines()[4:])+'\n'
  (d/'EXPEQ').write_text('\n'.join(lines)+'\n');(d/'chease_namelist').write_text(deck)
  cases[defect]=(d,None,None,None)
 results={};env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
 for name in ['CHEASE_GS_OBSERVER_DIR','CHEASE_NONLINEAR_OBSERVER_DIR','CHEASE_Q_OBSERVER_FILE']:env.pop(name,None)
 for name,(d,kappa,selector,expected)in cases.items():
  with(d/'chease_namelist').open('rb')as inp,(d/'parser.log').open('wb')as f:r=subprocess.run([str(build/'chease.x')],cwd=d,env=env,stdin=inp,stdout=f,stderr=subprocess.STDOUT,timeout=30)
  if selector is None:
   log=(d/'parser.log').read_text()
   assert r.returncode!=0,(name,'Malformed native input accepted')
   assert not(d/'parsed.bin').exists()
   assert 'Bad integer' in log if name=='missing_wall_header' else 'End of file' in log
   results[name]=dict(status='PASS: malformed input fails before PDE',returncode=r.returncode,log_sha256=sha(d/'parser.log'))
   continue
  assert r.returncode==0,(name,r.returncode)
  assert 'NATIVE INPUT ORACLE: NO PDE'in(d/'parser.log').read_text()
  assert not any((d/f).exists()for f in ['NOUT','NSAVE','EQDSK.OUT'])
  with FortranFile(d/'parsed.bin')as f:ints=f.read_ints('<i4');v=f.read_reals('<f8')
  assert ints.tolist()==[513,2,1,101,selector],(name,ints)
  t=np.arange(513)*2*np.pi/512
  assert np.max(abs(v[3:516]-(1+.1*np.cos(t))))<3e-16
  assert np.max(abs(v[516:1029]-kappa*.1*np.sin(t)))<2e-16
  assert np.max(abs(v[1029:1130]-np.linspace(0,1,101)))<2e-16
  assert np.all(v[1130:1231]==0)
  rf=np.full(101,1.5 if selector==4 else -.125)if expected is None else expected
  assert np.array_equal(v[1231:],rf),(name,'RFUN semantic readback')
  results[name]=dict(status='PASS',selector=int(ints[-1]),RFUN_min=float(rf.min()),RFUN_max=float(rf.max()),parsed_sha256=sha(d/'parsed.bin'),input_sha256={f:sha(d/f)for f in ['EXPEQ','chease_namelist']})
 receipt=dict(status='PASS',scope='Whole native INITIA/AUXVAL/BNDINP/IODISK; unconditional STOP in BNDINP immediately after IODISK returns. No PDE. Complete circle and shaped boundary, pressure, rho grid, selector and RFUN independently checked.',base_source_sha256=hashlib.sha256(original.encode()).hexdigest(),parser_source_sha256=sha(build/'chease.f'),binary_sha256=sha(build/'chease.x'),generator_sha256=sha(here/'generate_inputs.py'),script_sha256=sha(__file__),command=command,cases=results)
 (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':'PASS','cases':list(results),'receipt_sha256':sha(a.output/'receipt.json')}))
if __name__=='__main__':main()
