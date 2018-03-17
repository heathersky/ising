def T_anneal(T, step, num_steps, num_burnin):

    #implement annealing code here
    T_a = T
    T_0 = T_a+1
    m = -T_0/num_burnin
    T_step = m*step+T_0

    return float(max(T_step,T_a))

def B_anneal(B, step, num_steps, num_burnin):

    #implement annealing code here
    
    B_a = B
    B_0 = B_a+1
    s = -B_0/num_burnin
    B_step = s*step+B_0

    return float(max(B_step,B_a))