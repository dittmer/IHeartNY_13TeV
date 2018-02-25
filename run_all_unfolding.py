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

(options, args) = parser.parse_args()
argv = []

if options.closureTests and options.doSysClosure:
    path = [
        ## closure tests
        # Optional flags:
        ## --doSys to add systematic uncertainties
        ## --fullRange to extend pt range to [350,2000]
        ## --areaConstraint to add area constraint
        
        # Curvature regularization, [400,1200], no area constraint, ScanTau
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --doSys",
        # Curvature regularization, [400,1200], no area constraint, LCurve
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=curvature --doSys",
        # Derivative regularization, [400,1200], no area constraint, ScanTau
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=derivative --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=derivative --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=derivative --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=derivative --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=derivative --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=derivative --doSys",
        # Derivative regularization, [400,1200], no area constraint, LCurve
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=derivative --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=derivative --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=derivative --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=derivative --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=derivative --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=derivative --doSys",
        # Curvature regularization, [400,1200], area constraint, ScanTau
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --areaConstraint --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --areaConstraint --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=curvature --areaConstraint --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=curvature --areaConstraint --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --areaConstraint --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --areaConstraint --doSys",
        # Curvature regularization, [400,1200], area constraint, LCurve
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=curvature --areaConstraint --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=curvature --areaConstraint --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=curvature --areaConstraint --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=curvature --areaConstraint --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=curvature --areaConstraint --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=curvature --areaConstraint --doSys",
        # Curvature regularization, [350,2000], no area constraint, ScanTau
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --fullRange --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --fullRange --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=curvature --fullRange --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=curvature --fullRange --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --fullRange --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --fullRange --doSys",
        # Curvature regularization, [350,2000], no area constraint, LCurve
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=curvature --fullRange --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=curvature --fullRange --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=curvature --fullRange --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=curvature --fullRange --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=curvature --fullRange --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=curvature --fullRange --doSys",
        # Non-regularized unfolding
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=Tau0 --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=Tau0 --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=Tau0 --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=Tau0 --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=Tau0 --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=Tau0 --regMode=curvature --doSys",
        # Non-regularized unfolding, unfold full
        "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=nom --tauMode=Tau0 --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=nom --tauMode=Tau0 --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=up --tauMode=Tau0 --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=up --tauMode=Tau0 --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=dn --tauMode=Tau0 --regMode=curvature --doSys",
        "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=dn --tauMode=Tau0 --regMode=curvature --doSys",
    ]
elif options.closureTests:
    path = [
        # Curvature regularization, [400,1200], no area constraint, ScanTau
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=curvature",
        # Curvature regularization, [400,1200], no area constraint, LCurve
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=curvature",
        # Derivative regularization, [400,1200], no area constraint, ScanTau
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=derivative",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=derivative",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=derivative",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=derivative",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=derivative",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=derivative",
        # Derivative regularization, [400,1200], no area constraint, LCurve
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=derivative",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=derivative",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=derivative",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=derivative",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=derivative",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=derivative",
        # Curvature regularization, [400,1200], area constraint, ScanTau
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --areaConstraint",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --areaConstraint",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=curvature --areaConstraint",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=curvature --areaConstraint",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --areaConstraint",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --areaConstraint",
        # Curvature regularization, [400,1200], area constraint, LCurve
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=curvature --areaConstraint",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=curvature --areaConstraint",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=curvature --areaConstraint",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=curvature --areaConstraint",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=curvature --areaConstraint",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=curvature --areaConstraint",
        # Curvature regularization, [350,2000], no area constraint, ScanTau
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --fullRange",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=ScanTau --regMode=curvature --fullRange",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=ScanTau --regMode=curvature --fullRange",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=ScanTau --regMode=curvature --fullRange",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --fullRange",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=ScanTau --regMode=curvature --fullRange",
        # Curvature regularization, [350,2000], no area constraint, LCurve
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=LCurve --regMode=curvature --fullRange",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=LCurve --regMode=curvature --fullRange",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=LCurve --regMode=curvature --fullRange",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=LCurve --regMode=curvature --fullRange",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=LCurve --regMode=curvature --fullRange",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=LCurve --regMode=curvature --fullRange",
        # Non-regularized unfolding
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=nom --tauMode=Tau0 --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=nom --tauMode=Tau0 --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=up --tauMode=Tau0 --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=up --tauMode=Tau0 --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=muon --type=half --toy=dn --tauMode=Tau0 --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=half --toy=dn --tauMode=Tau0 --regMode=curvature",
        # Non-regularized unfolding, unfold full
        "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=nom --tauMode=Tau0 --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=nom --tauMode=Tau0 --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=up --tauMode=Tau0 --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=up --tauMode=Tau0 --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=muon --type=full --toy=dn --tauMode=Tau0 --regMode=curvature",
        "python closureTestTUnfold_v2.py --lepType=ele --type=full --toy=dn --tauMode=Tau0 --regMode=curvature",
    ]

else :
    path = [
        "python unfoldTopPt.py"
    ]

## run actual unfolding
for s in path :
    print s
    subprocess.call( [s, ""], shell=True )
    
    
