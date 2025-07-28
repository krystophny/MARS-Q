% retrieve D3D MDS+ data with Matlab from IRIS:
% for threshold I-coil current amplitude and phase for ELM suppression

% generate g-file from OMFIT
% - import EFIT (double click)
% - select "load existing g-file"
% - load from MDS+
% - EFIT Tree: CAKE01 or CAKE02
% - save project as ...
% - copy g-file from /fusion/projects/...

% generate pfile from OMFIT:
% - import kineticEFITtime
% - select "OMFITprofiles"
% - run Main
% - Export: to PFILES

% get DIII-D I-coil current:
% - use this Matlat script, or
% - run reviewplus
% - plot signals: time, pcil30-330, pciu30-330

shot = 179358;
time = 2560;    %time for threshold current
n    = '3';

mdsconnect('atlas.gat.com');
mdsopen('d3d',shot);

AU = mdsvalue(['\iun' n 'iamp']);  %max current in each coil
PU = mdsvalue(['\iun' n 'iphase']);
TU = mdsvalue(['dim_of(\iun' n 'iamp)']);

AL = mdsvalue(['\iln' n 'iamp']);
PL = mdsvalue(['\iln' n 'iphase']);
TL = mdsvalue(['dim_of(\iln' n 'iamp)']);

[T,IU] = min(abs(TU-time));
[T,IL] = min(abs(TL-time));

RES = [AU(IU) AL(IL) (AU(IU)+AL(IL))/2 PU(IU) PL(IL)]


