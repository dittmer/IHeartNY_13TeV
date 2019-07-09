## PARTON PT 
rm CrossSection_Parton_pt.root

python unfoldTopPt.py --norm
python unfoldTopPt.py

hadd CrossSection_Parton_pt.root CrossSection_Parton_pt_TMP.root CrossSection_Parton_pt_norm_TMP.root 
rm CrossSection_Parton_pt_TMP.root CrossSection_Parton_pt_norm_TMP.root 

#root -b -q 'DrawCrossSection.C("Parton","pt","p_{T}^{t} [GeV]","d#sigma/dp_{T}^{t} [pb/GeV]","1/#sigma d#sigma/dp_{T}^{t} [pb/GeV]",0.0001, 0.1, -2, 2, 0.0001, 0.012, -2, 2, false)' 
root -b -q 'DrawCrossSection.C("Parton","pt","p_{T}^{t} [GeV]","d#sigma/dp_{T}^{t} [pb/GeV]","1/#sigma d#sigma/dp_{T}^{t} [pb/GeV]",0.000008, 0.4, -1, 1, 0.000002, 0.04, -1, 1, true)'


## PARTICLE PT
rm CrossSection_Particle_pt.root

python unfoldTopPt.py --level="part" --norm
python unfoldTopPt.py --level="part"

hadd CrossSection_Particle_pt.root CrossSection_Particle_pt_TMP.root CrossSection_Particle_pt_norm_TMP.root 
rm CrossSection_Particle_pt_TMP.root CrossSection_Particle_pt_norm_TMP.root 

#root -b -q 'DrawCrossSection.C("Particle","pt","p_{T}^{t} [GeV]","d#sigma/dp_{T}^{t} [pb/GeV]","1/#sigma d#sigma/dp_{T}^{t} [pb/GeV]",0.0001, 0.025, -2, 2, 0.0001, 0.01, -2, 2, false)' 
root -b -q 'DrawCrossSection.C("Particle","pt","p_{T}^{t} [GeV]","d#sigma/dp_{T}^{t} [pb/GeV]","1/#sigma d#sigma/dp_{T}^{t} [pb/GeV]",0.00002, 0.06, -1, 1, 0.00001, 0.02, -1, 1, true)'

## PARTON RAPIDITY
rm CrossSection_Parton_y.root

python unfoldTopPt.py --toUnfold="y" --norm
python unfoldTopPt.py --toUnfold="y" 

hadd CrossSection_Parton_y.root CrossSection_Parton_y_TMP.root CrossSection_Parton_y_norm_TMP.root 
rm CrossSection_Parton_y_TMP.root CrossSection_Parton_y_norm_TMP.root 

#root -b -q 'DrawCrossSection.C("Parton","y","|y^{t}|","d#sigma/dy^{t} [pb]","1/#sigma d#sigma/dy^{t}",0, 10, -2, 2, 0, 1, -2, 2, false)' 
root -b -q 'DrawCrossSection.C("Parton","y","|y^{t}|","d#sigma/dy^{t} [pb]","1/#sigma d#sigma/dy^{t}",0, 8, -0.8, 0.8, 0, 1, -0.8, 0.8, false)' 

## PARTICLE RAPIDITY
rm CrossSection_Particle_y.root

python unfoldTopPt.py --level="part" --toUnfold="y" --norm
python unfoldTopPt.py --level="part" --toUnfold="y" 

hadd CrossSection_Particle_y.root CrossSection_Particle_y_TMP.root CrossSection_Particle_y_norm_TMP.root 
rm CrossSection_Particle_y_TMP.root CrossSection_Particle_y_norm_TMP.root 

#root -b -q 'DrawCrossSection.C("Particle","y","|y^{t}|","d#sigma/dy^{t} [pb]","1/#sigma d#sigma/dy^{t}",0, 5, -2, 2, 0, 1, -2, 2, false)' 
root -b -q 'DrawCrossSection.C("Particle","y","|y^{t}|","d#sigma/dy^{t} [pb]","1/#sigma d#sigma/dy^{t}",0, 3, -0.8, 0.8, 0, 1, -0.8, 0.8, false)' 

