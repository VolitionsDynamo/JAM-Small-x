def _hook1(params,order):
    sample = []
    for i in range(len(order)):
        if order[i][0] != 1: continue
        if order[i][1] != 'pdf': continue
        if 'g1'  in order[i][2]: continue
        if 'uv1' in order[i][2]: continue
        # if 'dv1' in order[i][2]: continue
        if 'db1' in order[i][2]: continue
        if 'ub1' in order[i][2]: continue
        if 's1' in order[i][2]: continue
        if 'sb1' in order[i][2]: continue
        sample.append(params[i])
    return sample


def _hook2(params,order):
    sample = []
    for i in range(len(order)):
        if order[i][0] != 1: continue
        if order[i][1] != 'pdf': continue
        if 'g1'  in order[i][2]: continue
        if 'uv1' in order[i][2]: continue
        # if 'dv1' in order[i][2]: continue
        if 'db1' in order[i][2]: continue
        if 'ub1' in order[i][2]: continue
        if 's1' in order[i][2]: continue
        if 'sb1' in order[i][2]: continue
        sample.append(params[i])
    return sample

def _hook3(params,order):
    sample = []
    for i in range(len(order)):
        if order[i][0] != 1: continue
        if order[i][1] != 'pdf': continue
        if 'g1'  in order[i][2]: continue
        if 'uv1' in order[i][2]: continue
        # if 'dv1' in order[i][2]: continue
        if 'db1' in order[i][2]: continue
        if 'ub1' in order[i][2]: continue
        if 's1' in order[i][2]: continue
        if 'sb1' in order[i][2]: continue
        sample.append(params[i])
    return sample

def _hook4(params,order):
    sample = []
    for i in range(len(order)):
        if order[i][0] != 1: continue
        if order[i][1] != 'pdf': continue
        if 'g1'  in order[i][2]: continue
        if 'uv1' in order[i][2]: continue
        if 'dv1' in order[i][2]: continue
        # if 'db1' in order[i][2]: continue
        if 'ub1' in order[i][2]: continue
        if 's1' in order[i][2]: continue
        if 'sb1' in order[i][2]: continue
        sample.append(params[i])
    return sample

def _hook5(params,order):
    sample = []
    for i in range(len(order)):
        if order[i][0] != 1: continue
        if order[i][1] != 'pdf': continue
        if 'g1'  in order[i][2]: continue
        if 'uv1' in order[i][2]: continue
        if 'dv1' in order[i][2]: continue
        # if 'db1' in order[i][2]: continue
        if 'ub1' in order[i][2]: continue
        if 's1' in order[i][2]: continue
        if 'sb1' in order[i][2]: continue
        sample.append(params[i])
    return sample

def _hook6(params,order):
    sample = []
    for i in range(len(order)):
        if order[i][0] != 1: continue
        if order[i][1] != 'pdf': continue
        if 'g1'  in order[i][2]: continue
        # if 'uv1' in order[i][2]: continue
        if 'dv1' in order[i][2]: continue
        if 'db1' in order[i][2]: continue
        if 'ub1' in order[i][2]: continue
        if 's1' in order[i][2]: continue
        if 'sb1' in order[i][2]: continue
        sample.append(params[i])
    return sample

def _hook7(params,order):
    sample=[]
    for i in range(len(order)):
        if order[i][0] != 1: continue
        if order[i][1] != 'ppdf': continue
        if 'g1'  in order[i][2]: continue
        if 'up1' in order[i][2]: continue
        if 'dp1' in order[i][2]: continue
        if 'sp1' in order[i][2]: continue
        sample.append(params[i])
    return sample

def _hook8(params,order):
    sample=[]
    for i in range(len(order)):
        if order[i][0] != 1: continue
        if order[i][1] != 'ppdf': continue
        # if 'g1'  in order[i][2]: continue
        if 'uv1' in order[i][2]: continue
        if 'dv1' in order[i][2]: continue
        if 'ub1' in order[i][2]: continue
        if 'db1' in order[i][2]: continue
        if 's1' in order[i][2]: continue
        if 'sb1' in order[i][2]: continue
        if 'sea1' in order[i][2]: continue
        sample.append(params[i])
    return sample

def _hook9(params,order):
    sample=[]
    for i in range(len(order)):
        if order[i][0] != 1: continue
        if order[i][1] != 'pdf': continue
        if 'g1'  in order[i][2]: continue
        if 'uv1' in order[i][2]: continue
        # if 'dv1' in order[i][2]: continue
        if 'db1' in order[i][2]: continue
        if 'ub1' in order[i][2]: continue
        if 's1' in order[i][2]: continue
        if 'sb1' in order[i][2]: continue

        # if 'g2'  in order[i][2]: continue
        # if 'uv2' in order[i][2]: continue
        # if 'dv2' in order[i][2]: continue
        # if 'db2' in order[i][2]: continue
        # if 'ub2' in order[i][2]: continue
        # if 's2' in order[i][2]: continue
        # if 'sb2' in order[i][2]: continue
        sample.append(params[i])
    return sample

def _hook10(params,order):
    sample=[]
    for i in range(len(order)):
        if order[i][0]!=1: continue
        #if order[i][1]!='pdf': continue
        #if 'g1'  in order[i][2]: continue
        #if 'uv1' in order[i][2]: continue
        #if 'dv1' in order[i][2]: continue
        #if 'db1' in order[i][2]: continue
        #if 'ub1' in order[i][2]: continue
        if 's1' in order[i][2]: continue
        #if 'sb1' in order[i][2]: continue
        sample.append(params[i])
    return sample

nc = {}
hooks = {}

## there is no need to cluster in first step
nc[1] = 1
hooks[1] = None

nc[2] = 1
hooks[2] = None

nc[3] = 1
hooks[3] = None

nc[4] = 1
hooks[4] = None

nc[5] = 1
hooks[5] = None

nc[6] = 1
hooks[6] = None

nc[7] = 2
hooks[7] = None

nc[8] = 1
hooks[8] = None

nc[9] = 1
hooks[9] = None

nc[10] = 1
hooks[10] = None
