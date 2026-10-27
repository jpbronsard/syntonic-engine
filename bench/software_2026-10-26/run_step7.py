"""Runs step 7 from the logs in data/ (gzip). The five step-7 scripts are byte-identical to the ones that ran on
27 Oct 2026 (their SHA-256 is listed in NOTE.md, section 12; the two pre-registered ones also print it at
the top of their results); they read the decoded logs from fixed paths of the
working session, so this runner unpacks data/ into a temporary folder and points them there.

data/ holds three sessions of the X500's SD logger, reduced to the columns the scripts use (time, armed, current,
throttle, raw accelerometer at 250 Hz, raw gyroscope at 1 kHz as four snapshots per record, and the engine's three
rate windows), copied as text without reformatting. No position or GPS column is included.
  sess6.csv     session 6: motor test with ballast, frame still, 0.5 to 1.2 A
  syntonic_00000000_000000_00005.csv   session 5: propulsion on a held frame, up to 0.44 throttle and 10 A
  crash7.csv    session 7: about 5 s airborne, then a throttle cut and a crash (cause outside the estimator)
"""
import contextlib, gzip, io, os, runpy, shutil, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp()
for f in ('sess6.csv', 'crash7.csv', 'syntonic_00000000_000000_00005.csv'):
    with gzip.open(os.path.join(HERE, 'data', f + '.gz'), 'rb') as a, open(os.path.join(tmp, f), 'wb') as b:
        shutil.copyfileobj(a, b)
sys.path.insert(0, HERE)

sys.argv = ['step7_axis_coupling.py', os.path.join(tmp, 'sess6.csv')]
import step7_axis_coupling as A                         # step 7a, pre-registered: attitude, roll and pitch
A.main()

buf = io.StringIO()                                     # exploratory: is the V5.1 window the best symmetric one?
with contextlib.redirect_stdout(buf):
    runpy.run_path(os.path.join(HERE, 'step7_exploratory_symmetric.py'), run_name='__main__')
print(buf.getvalue(), end='')
open(os.path.join(HERE, 'results', 'step7_exploratory_symmetric.txt'), 'w').write(buf.getvalue())

import step7c_rate3_th44 as C                           # step 7c, pre-registered: rates, roll, pitch, yaw, TH44
C.UP = tmp + os.sep
C.main()
runpy.run_path(os.path.join(HERE, 'step7c_A1_exploratory.py'), run_name='__main__')    # amendment A1
runpy.run_path(os.path.join(HERE, 'step7c_A1_sensitivity.py'), run_name='__main__')    # its sensitivity
shutil.rmtree(tmp)
