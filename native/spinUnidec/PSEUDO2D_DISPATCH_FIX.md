Pseudo2D restrained-fit dispatch fix
====================================

1. decon_parallel.cpp now mirrors decon.cpp: when dim==2 and pseudo2DFit==1,
   dispatch to Protocol2PFit rather than Protocol2D.
2. Protocol2PFit no longer rejects pseudo2D NMRPipe files whose header reports
   dim==3 solely because a singleton Z storage axis is retained. It validates
   the two active dimord axes instead.
3. Existing Protocol2D/Protocol3P paths are unchanged.
