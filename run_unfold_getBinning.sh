root -b -q 'unfold_getBinning.C("mu","pt",false)' 
root -b -q 'unfold_getBinning.C("el","pt",false)' 
root -b -q 'unfold_getBinning.C("mu","pt",true)' 
root -b -q 'unfold_getBinning.C("el","pt",true)' 

root -b -q 'unfold_getBinning.C("mu","y",false)' 
root -b -q 'unfold_getBinning.C("el","y",false)' 
root -b -q 'unfold_getBinning.C("mu","y",true)' 
root -b -q 'unfold_getBinning.C("el","y",true)' 

root -b -q drawCombineMatrix.C