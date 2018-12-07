import subprocess
import sys

from optparse import OptionParser
parser = OptionParser()

parser.add_option('--closureTests', metavar='F', action='store_true',
                  default=False,
                  dest='closureTests',
                  help='Run closure test')

parser.add_option('--doSysClosure', metavar='F', action='store_true',
                  default=False,
                  dest='doSysClosure',
                  help='Use systematic uncertainties in closure test')

parser.add_option('--toUnfold', metavar='F', type='string', action='store',
                  default='pt',
                  dest='toUnfold',
                  help='Distribution to unfold (pt or y)')

parser.add_option('--level', metavar='F', type='string', action='store',
                  default='gen',
                  dest='level',
                  help='Level to unfold (gen or part)')

(options, args) = parser.parse_args()
argv = []

if options.closureTests :
    path = [
        ## closure tests
        # Optional flags:
        ## --doSys to add systematic uncertainties
        ## --fullRange to extend pt range to [350,2000]
        
        # Curvature regularization, [400,1200], ScanTau
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        # Curvature regularization, [400,1200], LCurve
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        # Derivative regularization, [400,1200], ScanTau
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        # Derivative regularization, [400,1200], LCurve
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=derivative --toUnfold="+options.toUnfold+" --level="+options.level,
        # Non-regularized unfolding
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        # Non-regularized unfolding, unfold full
        "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=nom --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=nom --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=up --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=up --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=dn --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
        "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=dn --tauMode=Tau0 --regMode=curvature --toUnfold="+options.toUnfold+" --level="+options.level,
    ]
    
    if options.toUnfold == "pt":
        path += [
            # Curvature regularization, [350,2000], ScanTau
            "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            # Curvature regularization, [350,2000], LCurve
            "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            # Derivative regularization, [350,2000], ScanTau
            "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            # Derivative regularization, [350,2000], LCurve
            "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=derivative --fullRange --toUnfold=pt --level="+options.level,
            # Non-regularized unfolding
            "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            #"python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            # Non-regularized unfolding, unfold full
            "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=nom --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=nom --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=up --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=up --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=dn --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
            "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=dn --tauMode=Tau0 --regMode=curvature --fullRange --toUnfold=pt --level="+options.level,
        ]
        

else :
    path = [
        "python unfoldTopPt.py --toUnfold="+options.toUnfold+" --level="+options.level
    ]

## run actual unfolding
for s in path :
    if options.closureTests and options.doSysClosure :
        s += " --doSys"
    print s
    subprocess.call( [s, ""], shell=True )
    
    
