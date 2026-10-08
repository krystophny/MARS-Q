"""Hermite chart interpolation alone creates inverse-q axis conditioning."""
import hashlib
import itertools
import json
from pathlib import Path
import argparse

import mpmath as mp
import numpy as np



def local_polynomial(state, point):
    # Exact local bicubic basis in the literal circular chart.
    sigma = np.linalg.norm(point-np.array([1.,0.]))/.1
    theta = np.arctan2(point[1],point[0]-1.)
    s, t = state['radial_knots'], state['angular_knots']
    i = np.searchsorted(s,sigma)-1;j = np.searchsorted(t,theta)-1
    a = state['coefficients'].reshape(len(s),len(t)-1,4)
    sm,tm,hs,ht = map(mp.mpf,map(float,[s[i],t[j],s[i+1]-s[i],t[j+1]-t[j]]))
    c = [[[mp.mpf(float(x))for x in a[i+ii,(j+jj)%(len(t)-1)]]
         for jj in range(2)]for ii in range(2)]
    def basis(x,h):return [2*x**3-3*x*x+1,h*(x**3-2*x*x+x),-2*x**3+3*x*x,h*(x**3-x*x)]
    def value(R,Z):
        sv=basis((mp.sqrt((R-1)**2+Z*Z)/mp.mpf('.1')-sm)/hs,hs)
        tv=basis((mp.atan2(Z,R-1)-tm)/ht,ht)
        v=mp.mpf(0)
        for ii in range(2):
            for jj in range(2):
                cc=c[ii][jj]
                v+=cc[0]*sv[2*ii]*tv[2*jj]+cc[1]*sv[2*ii+1]*tv[2*jj]+cc[2]*sv[2*ii]*tv[2*jj+1]+cc[3]*sv[2*ii+1]*tv[2*jj+1]
        return v
    return value,dict(radial_element=int(i),angular_element=int(j),
        sigma=sigma,theta=theta,sigma_interval=[float(s[i]),float(s[i+1])],
        theta_interval=[float(t[j]),float(t[j+1])])


def coarea_jets(R,H,T,U):
    eig,vec=np.linalg.eigh(H);assert eig.min()>0
    A=(vec/np.sqrt(eig))@vec.T
    theta=np.arange(2048)*2*np.pi/2048
    n=np.c_[np.cos(theta),np.sin(theta)]
    e=n@A.T
    c=np.einsum('ijk,ni,nj,nk->n',T,e,e,e)/6
    d=np.einsum('ijkl,ni,nj,nk,nl->n',U,e,e,e,e)/24
    terms=np.stack([e[:,0]**2/R**3,4*c*e[:,0]/R**2,
                    12*c*c/R,-4*d/R])
    K0=2*np.pi/(R*np.sqrt(np.linalg.det(H)))
    contributions=4*np.pi/np.sqrt(np.linalg.det(H))*terms.mean(axis=1)
    Kdelta=contributions.sum()
    D2delta=-2*Kdelta/K0**3
    return dict(K_axis=K0,D_axis=1/K0,K_derivative_flux=Kdelta,
                K_derivative_terms=contributions.tolist(),
                terms=['R_weight','cubic_weight','cubic_squared','quartic'],
                D2_derivative_flux=D2delta,
                constant_q_FFprime=(2*np.pi*1.5)**2*D2delta/2)



def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path)
    a=p.parse_args();assert not a.output.exists()
    mp.mp.dps=55;C=10.;Ra=1.00125;b=.1;rows=[]
    for nt in [24,48,96,192]:
        s=np.linspace(0,1,25);t=np.linspace(-np.pi/nt,2*np.pi-np.pi/nt,nt+1)
        ss,tt=np.meshgrid(s,t[:-1],indexing='ij');x=b*ss*np.cos(tt)+1-Ra;y=b*ss*np.sin(tt)
        values=np.stack([-.1+C*(x*x+y*y),
            2*C*(b*b*ss+b*(1-Ra)*np.cos(tt)),
            -2*C*b*ss*(1-Ra)*np.sin(tt),
            -2*C*b*(1-Ra)*np.sin(tt)],axis=-1)
        state=dict(radial_knots=s,angular_knots=t,coefficients=values.ravel())
        f,location=local_polynomial(state,np.array([Ra,0.]))
        R=mp.findroot(lambda R:mp.diff(f,(R,mp.mpf(0)),(1,0)),mp.mpf(str(Ra)))
        Z=mp.mpf(0);jets={}
        for order in [2,3,4]:
            v=np.zeros((2,)*order)
            for indices in itertools.product(range(2),repeat=order):
                v[indices]=float(mp.diff(f,(R,Z),(indices.count(0),indices.count(1))))
            jets[order]=v
        observed=coarea_jets(float(R),jets[2],jets[3],jets[4])
        exact=coarea_jets(Ra,2*C*np.eye(2),np.zeros((2,)*3),np.zeros((2,)*4))
        rows.append(dict(NT=nt,interpolated_axis_R=float(R),exact_axis_R=Ra,
            axis_position_error=float(R)-Ra,chart=location,
            Hessian_max_absolute_error=float(np.max(abs(jets[2]-2*C*np.eye(2)))),
            Hermite_axis_FFprime=observed['constant_q_FFprime'],
            exact_physical_quadratic_FFprime=exact['constant_q_FFprime'],
            K_derivative_terms=observed['K_derivative_terms']))
    assert rows[0]['Hermite_axis_FFprime']>1e6
    d=dict(scope='Held interpolation of an exact shifted physical quadratic, using all analytically consistent value/radial/angular/mixed jets on the native chart and grids. This is not a boundary-value equilibrium or a PDE solve. Axis, C and radius are fixed analytic parameters; no fit to native source or field.',
        controls=rows,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        Taylor_reader_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    a.output.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps(rows))


if __name__=='__main__':main()
