"""Seconds-only boundary fixture; no model solve."""
import numpy as np
from strict_tenure import tenure_logit_kernel
b=np.array([-1.,0.,1.]); V=np.zeros((3,2,1,1,1)); heq=np.array([[0.,.5]])
h=np.zeros((1,2)); dp=np.zeros((1,2,1,1)); bm=np.full((1,2,1,1),-99.); bd=np.zeros((1,1,2,2),dtype=np.bool_); grant=np.zeros((1,2,1,1))
_,choice,prob=tenure_logit_kernel(V,b,heq,h,dp,bm,bd,grant,.1,V)
# owner sale balances: -0.5, 0.5, 1.5.  Negative equity cannot select renter.
assert choice[0,1,0,0,0] != 0 and prob[0,1,0,0,0,0] == 0
assert prob[1,1,0,0,0,0] > 0 and prob[2,1,0,0,0,0] > 0
print('strict sale fixture passed')
