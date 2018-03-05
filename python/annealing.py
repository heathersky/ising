def T_anneal(T, step, num_steps, num_burnin):

    #implement annealing code here
    T_a = float(T)
    T_0 = 100
    m = -T_0/float(num_burnin)
    T_step = m*float(step)+T_0

    return float(max(T_step,T_a))

def B_anneal(B, step, num_steps, num_burnin):

    #implement annealing code here
    
    B_a = B

    return float(B_a)