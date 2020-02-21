#!/usr/bin/env python
# coding: utf-8

# Functions used in the modified BRD 2D-ILT procedure

# In[1]:


pwd


# In[2]:


#References
#Venkataramanan et al. - 2002 (DOI : 10.1109/78.995059)
#Mitchell et al. - 2012 (DOI : 10.1016/j.pnmrs.2011.07.002)
def nnls_reg(val):
    r'''A function that initializes the storage arrays and performs nnls using smoothing parameter :val:
        Calls function :A_prime: which generates the smoothing matrix for the fit and :b_prime: which is
        the lexicographically ordered data, extended by zeros to match the dimensions of the fit matrix.
    Parameters
    ----------
    val : float
        The smoothing parameter lambda (= sqrt(alpha))
    '''
    r_norm = empty([1])
    x,r_norm[0] = nnls(A_prime(val,dimension),b_prime)
    return x

def chi(x_vec,val):
    r'''The function that minimizes the c_r vector derived from Kuhn-Tucker conditions.
        Eq[31] in Venkataramanan et al. - 2002.
    Parameters
    ----------
    x_vec : array
        The c_r vector
    val : float
        The smoothing parameter lambda, which is converted to alpha (i.e., lambda**2)
        in the function.
    '''
    return 0.5*dot(x_vec.T,dot(dd_chi(G(x_vec),val**2),x_vec)) - dot(x_vec.T,m_vec)

def d_chi(x_vec,val):
    r'''The derivative of :chi:, which serves as the input function to the Newton Minimization
        procedure -- i.e., we find the zeros of this function which correspond to maxima or minima
        of :chi:
        Eq[32] in Venkataramanan et al. - 2002.
    Parameters
    ----------
    x_vec : array
        The c_r vector
    val : float
        The smoothing parameter lambda, which is converted to alpha (i.e., lambda**2)
        in the function.
    '''
    return dot(dd_chi(G(x_vec),val**2),x_vec) - m_vec

def dd_chi(G,val):
    r'''The second derivative of :chi:, which serves as the derivative of the input function to the
        Newton Minimization procedure -- i.e., we find the zeros of this function which correspond to
        maxima or minima of :chi:
        Eq[33] in Venkataramanan et al. - 2002.
    Parameters
    ----------
    G : matrix
        The diagonal matrix, see :G: 
    val : float
        The smoothing parameter lambda, which is converted to alpha (i.e., lambda**2)
        in the function.
    '''
    return G + (val**2)*eye(shape(G)[0])

def G(x_vec):
    r'''The symmetric matrix used in the BRD algorithm.   
        Eq[30] in Venkataramanan et al. - 2002.
    Parameters
    ----------
    x_vec : array
        The c_r vector
    '''
    return dot(K0,dot(square_heaviside(x_vec),K0.T))

def H(product):
    r'''A simple heaviside function, used in :G:
        See eq[30] in Venkataramanan et al. - 2002.
    Parameters
    ----------
    product : float
        The product of dotting a row element of the tensor product kernel with the
        c_r vector, used in :G:.
        See eq[30] in Venkataramanan et al. - 2002.
    '''
    if product <= 0:
        return 0
    if product > 0:
        return 1

def square_heaviside(x_vec):
    r'''The diagonal matrix used in :G:
        See eq[30] in Venkataramanan et al. - 2002.
    Parameters
    ----------
    x_vec : array
        The c_r vector
    '''
    diag_heavi = []
    for q in range(shape(K0.T)[0]):
        pull_val = dot(K0.T[q,:],x_vec)
        temp = pull_val[0]
        temp = H(temp)
        diag_heavi.append(temp)
    diag_heavi = array(diag_heavi)
    square_heavi = diag_heavi*eye(shape(diag_heavi)[0])
    return square_heavi

def newton_min(input_vec,val):
    r'''The Newton-Raphson minimization algorithm, provided by scipy.optimize.newton,
        which needs to be redefined to enable the use of vectors and matrices in the
        minimization procedure. Performs a single minimization to optimize the c_r
        vector.
    Parameters
    ----------
    input_vec : array
        The c_r vector generated from the NNLS output.
        Eq[26] in Venkataramanan et al. - 2002.
    val : float
        The smoothing parameter lambda, which is converted to alpha (i.e., lambda**2)
        in the function.
    '''
    fder = dd_chi(G(input_vec),val)
    fval = d_chi(input_vec,val)
    newton_step = dot(linalg.inv(fder),fval)
    update_vec = input_vec + newton_step
    return update_vec

def optimize_alpha(input_vec,val):
    r'''Optimizes the smoothing parameter using eq[48] in Venkataramanan et al. - 2002.
    Parameters
    ----------
    input_vec : array
        The optimized c_r vector generated from Newton Minimization procedure.
        Eq[26] in Venkataramanan et al. - 2002.
    val : float
        The smoothing parameter lambda, which is converted to alpha (i.e., lambda**2)
        in the function.
    '''
    alpha_converged = False
    fac = sqrt(choose_s1*choose_s2)
    T = linalg.inv(dd_chi(G(input_vec),val**2))
    dot_prod = dot(input_vec.T,dot(T,input_vec))
    ans = dot_prod*fac
    ans = ans/linalg.norm(input_vec)
    ans = ans/(dot_prod)
    tol = 1e-3
    if abs(ans-val**2) <= tol:
        print("ALPHA HAS CONVERGED")
        alpha_converged = True
        return ans,alpha_converged
    return ans,alpha_converged

def mod_BRD(guess,maxiter=20):
    r'''The modified BRD method presented in Venkataramanan et al. - 2002
        for optimizing the smoothing parameter using to regualrize the NNLS fit.
    Parameters
    ----------
    guess : float
        The initial smoothing parameter as lambda - i.e., sqrt(alpha). Algorithm
        should converge within a few steps irrespective of this guess.
    maxiter : int
        The maximum number of iterations for the algorithm, set to 20 by default.
        Should not need to be more than this.
     '''
    smoothing_param = guess
    alpha_converged = False
    for iter in range(maxiter):
        print("*** *** ITERATION NO.",iter,"*** ***")
        print("*** CURRENT LAMBDA",smoothing_param," *** ")
        x_norm = empty([1])
        r_norm = empty([1])
        soln,r_norm = nnls(A_prime(smoothing_param,dimension),b_prime)
        f_vec = soln[:,newaxis]
        alpha = smoothing_param**2
        c_vec = dot(K0,f_vec) - m_vec
        c_vec /= -alpha
        new_c = newton_min(c_vec,smoothing_param)
        new_alpha,alpha_converged = optimize_alpha(new_c,smoothing_param)
        new_lambda = sqrt(new_alpha[0,0])
        if alpha_converged:
            print("*** OPTIMIZED LAMBDA",new_lambda," *** ")
            break  
        if not alpha_converged:
            print("*** UPDATED LAMBDA",new_lambda," *** ")
            smoothing_param = new_lambda
        if iter == maxiter-1:
            print("DID NOT CONVERGE.")
    return new_lambda


# Initializing dataset

# In[3]:


from pyspecdata import *
from scipy.optimize import leastsq,minimize,basinhopping,nnls
fl = figlist_var()
for date,id_string in [
        ('200115','CPMG_DNP_1')
        ]:
    filename = date+'_'+id_string+'.h5'
    nodename = 'signal'
    s = nddata_hdf5(filename+'/'+nodename,
            directory = getDATADIR(
                exp_type = 'test_equip'))
            #{{{ pulling acq params
    SW_kHz = s.get_prop('acq_params')['SW_kHz']
    nPoints = s.get_prop('acq_params')['nPoints']
    nEchoes = s.get_prop('acq_params')['nEchoes']
    nPhaseSteps = s.get_prop('acq_params')['nPhaseSteps']
    nScans = s.get_prop('acq_params')['nScans']
    p90_s = s.get_prop('acq_params')['p90_us']*1e-6
    deadtime_s = s.get_prop('acq_params')['deadtime_us']*1e-6
    deblank_s = s.get_prop('acq_params')['deblank_us']*1e-6
    marker_s = s.get_prop('acq_params')['marker_us']*1e-6
    tau1_s = s.get_prop('acq_params')['tau1_us']*1e-6
    pad_start_s = s.get_prop('acq_params')['pad_start_us']*1e-6
    pad_end_s = s.get_prop('acq_params')['pad_end_us']*1e-6
    #}}}
    orig_t = s.getaxis('t')
    acq_time_s = orig_t[nPoints]
    s.set_units('t','s')
    twice_tau = deblank_s + 2*p90_s + deadtime_s + pad_start_s + acq_time_s + pad_end_s + marker_s
    t2_axis = linspace(0,acq_time_s,nPoints)
    tE_axis = r_[1:nEchoes+1]*twice_tau
    s.setaxis('t',None)
    #s.setaxis('nScans',r_[0:nScans])
    s.chunk('t',['ph1','tE','t2'],[nPhaseSteps,nEchoes,-1])
    s.setaxis('ph1',r_[0.,2.]/4)
    s.setaxis('tE',tE_axis)
    s.setaxis('t2',t2_axis)
    fl.next(id_string+'raw data - chunking')
    fl.image(s)
    s.ft('t2', shift=True)
    fl.next(id_string+'raw data - chunking ft')
    fl.image(s)
    s.ft(['ph1'])
    s = s['ph1',1].C
    fl.next(id_string+' image plot coherence-- ft ')
    fl.image(s)
    s.ift('t2')
    fl.next(id_string+' image plot coherence ')
    fl.image(s)
    #s.mean('nScans',return_error=False)
    s.reorder('t2',first=True)
    t2_max = zeros_like(s.getaxis('power'))
    for x in range(len(s.getaxis('power'))):
        t2_max[x] = abs(s['power',x]['tE',0]).argmax('t2',raw_index=True).data
    s.setaxis('t2',lambda t: t -s.getaxis('t2')[int(t2_max.mean())])
    s.rename('tE','nEchoes').setaxis('nEchoes',r_[1:nEchoes+1])
    s.reorder('nEchoes',first=True)
    s.ft('t2')
    # as of right now, this plot doesn't show anything meaningful
    fl.next('check center')
    fl.image(s)
    f_axis = s.fromaxis('t2')
    def costfun(p):
        zeroorder_rad,firstorder = p
        phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
        phshift *= exp(-1j*2*pi*zeroorder_rad)
        corr_test = phshift * s
        return (abs(corr_test.data.imag)**2)[:].sum()
    iteration = 0
    def print_fun(x, f, accepted):
        global iteration
        iteration += 1
        #print (iteration, x, f, int(accepted))
        return
    sol = basinhopping(costfun, r_[0.,0.],
            minimizer_kwargs={"method":'L-BFGS-B'},
            callback=print_fun,
            stepsize=100.,
            niter=100,
            T=1000.
            )
    zeroorder_rad, firstorder = sol.x
    phshift = exp(-1j*2*pi*f_axis*(firstorder*1e-6))
    phshift *= exp(-1j*2*pi*zeroorder_rad)
    s *= phshift
    #print("RELATIVE PHASE SHIFT WAS {:0.1f}\us and {:0.1f}$^\circ$".format(
    #        firstorder,angle(zeroorder_rad)/pi*180))
    #if s['nEchoes':0]['t2':0]['power',-4].item().real > 0:
        #print s['nEchoes':0]['t2':0]['power',-4].item().real
        #print "Sign correction"
        #s *= -1
        #print s['nEchoes':0]['t2':0]['power',-4].item().real
    #print ndshape(s)
    s.reorder('power',first=True)
    fl.next('after phased - real ft')
    fl.image(s.real)
    fl.next('after phased - imag ft')
    fl.image(s.imag)
    power_axis_dBm = array(s.get_prop('meter_powers'))
    power_axis_W = zeros_like(power_axis_dBm)
    power_axis_W[:] = (1e-2*10**((power_axis_dBm[:]+10.)*1e-1))


# In[4]:


power_axis_W = r_[0,power_axis_W]


# In[5]:


enh = s['t2':(-150,200)].C.sum('t2')
print(ndshape(enh['nEchoes',1]))


# In[6]:


enh_ = enh['nEchoes',0]*-1


# In[7]:


figure()
plot(enh_,'.')


# In[8]:


print(ndshape(s))
image(s)


# In[9]:


small_s = s.C
image(small_s['t2':(-100,150)])
d = small_s['power',:-3]['t2':(-100,150)].C.sum('t2')*-1


# In[10]:


image(d)
print(ndshape(d))


# In[13]:


# Get dataset of the indirect dimensions
# Create grid of parameter values
# Carefully define kernels and solve


# In[18]:


print(ndshape(d))
print(d.getaxis('power'))
print(d.getaxis('nEchoes'))


# In[20]:


d.setaxis('nEchoes',tE_axis)
print(ndshape(d))
print(d.getaxis('power')) # units of Watts
print(d.getaxis('nEchoes')) # units of sec


# In[65]:


# Create the power axis
pp5 = 0.3
pp5_0 = linspace(1e-3,pp5,25)
pp5_1 = linspace(pp5,1.1,25)
pp5_ax = ones(50)
pp5_ax[:25] = pp5_0[:]
pp5_ax[25:] = pp5_1[:]
pp5_ax = nddata(pp5_ax,'pp5')


# In[66]:


# Create the T2 axis
Ny = 50
T2_ax = nddata(logspace(-2.5,1,Ny),'T2')


# In[67]:


x = d.C.nnls(('power','nEchoes'),
            (pp5_ax,T2_ax),
            (lambda x1,x2: -x1/(x1+x2),
            lambda y1,y2: exp(-y1/y2)),
            l='BRD')


# In[68]:


x.setaxis('T2',log10(T2_ax.data)).set_units('T2',None)
x.set_units('pp5','W')


# In[69]:


image(x)


# In[70]:


x.contour(labels=False)


# In[71]:


fit = x.C.get_prop('K1').dot(x.C.dot(x.C.get_prop('K2')))


# In[73]:


residual  = d - fit


# In[77]:


figure(figsize=(6,15));suptitle('DATASET: %s_%s'%(date,id_string))
subplot(311);subplot(311).set_title('DATA\n $ m $')
image(d)
#maxval=abs(d.data).max()
#clim(-maxval,maxval)

subplots_adjust(hspace=1.0)
subplot(312);subplot(312).set_title('FIT\n $K_{1}$ $f$ $K_{2}$')
image(fit)
#clim(-maxval,maxval)

subplots_adjust(hspace=1.0)
subplot(313);subplot(313).set_title('RESIDUAL\n $ m $ - $K_{1}$ $f$ $K_{2}$')
image(residual)
#clim(-maxval,maxval)
savefig('20200219_T2DNP_resid.pdf',transparent=True,bbox_inches='tight',pad_inches=0)


# In[80]:


figure(figsize=(6,15));suptitle('DATASET: %s_%s'%(date,id_string))
subplot(311);subplot(311).set_title('DATA\n $ m $')
image(abs(d))
#maxval=abs(d.data).max()
#clim(-maxval,maxval)

subplots_adjust(hspace=1.0)
subplot(312);subplot(312).set_title('FIT\n $K_{1}$ $f$ $K_{2}$')
image(abs(fit))
#clim(-maxval,maxval)

subplots_adjust(hspace=1.0)
subplot(313);subplot(313).set_title('RESIDUAL\n $ m $ - $K_{1}$ $f$ $K_{2}$')
image(abs(residual))
#clim(-maxval,maxval)
savefig('20200219_T2DNP_resid_abs.pdf',transparent=True,bbox_inches='tight',pad_inches=0)


# In[11]:


data_nd = d.C
data_nd.setaxis('nEchoes',tE_axis)
data_nd_mesh = data_nd.C
data_nd_mesh *= -1
data_nd_mesh.meshplot(cmap=cm.viridis)
show()
data = data_nd.data


# In[ ]:


pp5 = 0.5
pp5_0 = linspace(1e-3,pp5,25)
pp5_1 = linspace(pp5,0.9,25)
pp5_ax = ones(50)
pp5_ax[:25] = pp5_0[:]
pp5_ax[25:] = pp5_1[:]


# In[ ]:


print("Constructing kernels...")
Nx = 50
Ny = 50
log_Nx_ax = pp5_ax # p1/2
#log_Ny_ax = linspace(log10(3e-3),log10(3.0),Ny) # T2
log_Ny_ax = linspace(log10(1e-1),log10(1.0),Ny) # T2
tau1 = data_nd.getaxis('power')
tau2 = data_nd.getaxis('nEchoes')
N1_4d = reshape(tau1,(shape(tau1)[0],1,1,1))
N2_4d = reshape(tau2,(1,shape(tau2)[0],1,1))
Nx_4d = reshape(log_Nx_ax,(1,1,shape(log_Nx_ax)[0],1))
Ny_4d = reshape(10**log_Ny_ax,(1,1,1,shape(log_Ny_ax)[0]))
print("Shape of parameter of interest (x) axis",shape(log_Nx_ax),shape(Nx_4d))
print("Shape of parameter of interest (y) axis",shape(log_Ny_ax),shape(Ny_4d))
print("Shape of indirect dimension (tau1) axis",shape(tau1),shape(N1_4d))
print("Shape of indirect dimension (tau2) axis",shape(tau2),shape(N2_4d))


# In[ ]:


#k1 = (1.-2*exp(-N1_4d/Nx_4d))
k1 = -N1_4d/(N1_4d+Nx_4d)
k2 = exp(-N2_4d/Ny_4d)
print("Shape of K1 (relates tau1 and x)",shape(k1))
print("Shape of K2 (relates tau2 and y)",shape(k2))
k1_sqz = squeeze(k1)
k2_sqz = squeeze(k2)
U1,S1_row,V1 = np.linalg.svd(k1_sqz,full_matrices=False)
print("SVD of K1",map(lambda x: x.shape, (U1, S1_row, V1)))
U2,S2_row,V2 = np.linalg.svd(k2_sqz,full_matrices=False)
print("SVD of K2",map(lambda x: x.shape, (U2, S2_row, V2)))


# In[ ]:


figure()
semilogy(S1_row,'o-',alpha=0.4)
semilogy(S2_row,'o-',alpha=0.4)


# In[ ]:


print("")
print("*** BEGINNING COMPRESSION ***")
print("")
data_max = amax(data_nd.data)
print("Maximum in the data",data_max)
choose_s1 = 8
choose_s2 = 8


# In[ ]:


print("Uncompressed singular row vector for K1",S1_row.shape)
S1_row = S1_row[0:choose_s1]
print("Compressed singular value row vector for K1",S1_row.shape)
V1 = V1[0:choose_s1,:]
U1 = U1[:,0:choose_s1]
print("Compressed V matrix for K1",V1.shape)
print("Comrpessed U matrix for K1",U1.shape)

print("Uncompressed singular row vector for K2",S2_row.shape)
S2_row = S2_row[0:choose_s2]
print("Compressed singular value row vector for K2",S2_row.shape)
V2 = V2[0:choose_s2,:]
U2 = U2[:,0:choose_s2]
print("Compressed V matrix for K2",V2.shape)
print("Compressed U matrix for K2",U2.shape)

I_S1 = eye(S1_row.shape[0])
S1 = S1_row*I_S1
print("Non-zero singular value matrix for K1",S1.shape)

I_S2 = eye(S2_row.shape[0])
S2 = S2_row*I_S2
print("Non-zero singular value matrix for K2",S2.shape)


# In[ ]:


data_proj = U1.dot(U1.T.dot(data.dot(U2.dot(U2.T))))

for tau1_index in range(shape(data_proj)[0]):
    title('projected data')
    plot(data_proj[tau1_index,:])
show()


# In[ ]:


nd_proj = nddata(data_proj,['N1','N2'])
nd_proj.name('Projected data')
nd_proj.setaxis('N1',data_nd.getaxis('power')).rename('N1',r'$\tau_{1}$')
nd_proj.setaxis('N2',data_nd.getaxis('nEchoes')).rename('N2',r'$\tau_{2}$')
nd_proj.meshplot(cmap=cm.viridis)


# In[ ]:


print("Projected data dimensions:",shape(data_proj))
data_compr = U1.T.dot(data.dot(U2))
print("Compressed data dimensioins:",shape(data_compr))


# In[ ]:


comp = data_compr
comp = reshape(comp,(shape(data_compr))[0]*(shape(data_compr))[1])


# In[ ]:


figure()
for x in range((shape(data_compr))[1]):
    plot(data_compr[:,x],'-.')
ylabel('Compressed data')
xlabel('Index')
show()


# In[ ]:


figure()
plot(comp,'-.')
ylabel('Compressed data')
xlabel('Index')
show()


# In[ ]:


K1 = S1.dot(V1)
K2 = S2.dot(V2)
print("Compressed K1",shape(K1))
print("Compressed K2",shape(K2))


# In[ ]:


K1 = reshape(K1, (shape(K1)[0],1,shape(K1)[1],1))
K2 = reshape(K2, (1,shape(K2)[0],1,shape(K2)[1]))
K0 = K1*K2
K0 = reshape(K0, (shape(K1)[0]*shape(K2)[1],shape(K1)[2]*shape(K2)[3]))
print("Compressed tensor kernel",shape(K0))
print("* Should be (",shape(S1)[0],"*",shape(S2)[0],") x (",shape(Nx_4d)[2],"*",shape(Ny_4d)[3],")")
#END COMPRESSION


# In[ ]:


datac_lex = []
for m in range(shape(data_compr)[0]):
    for l in range(shape(data_compr)[1]):
        temp = data_compr[m][l]
        datac_lex.append(temp)
print("Dimension of lexicographically ordered data:",shape(datac_lex)[0])
print("Should match first dimension of compressed tensor kernel K0 which is",shape(K0)[0])


# In[ ]:


x, rnorm = nnls(K0,datac_lex)
solution = reshape(x,(Nx,Ny))
figure()
title('Estimate, no regularization')
image(solution)
show()
   


# In[ ]:


datac_lex = array(datac_lex)
datac_lex = datac_lex[:,newaxis]
print("Lexicographically orderd data:",shape(datac_lex))


# In[ ]:


dimension = K0.shape[1]
def A_prime(val,dimension):
    A_prime = r_[K0, val*eye(dimension)]
    return A_prime

b_prime = r_[datac_lex,zeros((dimension,1))]
b_prime = b_prime.squeeze()
print("Shape of b vector",shape(b_prime))
m_vec = datac_lex


# In[ ]:


if S_curve:
    if gen_S_curve_data:
        lambda_range = logspace(log10(8e-4),log10(2e4),3)
        rnorm_list = empty_like(lambda_range)
        smoothing_list = empty_like(lambda_range)
        alt_norm_list = empty_like(lambda_range)
        for index,lambda_val in enumerate(lambda_range):
            print "index",index
            soln,temp_rn = nnls(A_prime(lambda_val,dimension),b_prime)
            rnorm_list[index] = temp_rn
            f_vec = soln[:,newaxis]
            alpha = lambda_val**2
            c_vec = dot(K0,f_vec) - m_vec
            c_vec /= -alpha
            alt_temp = linalg.norm(c_vec)*alpha
            alt_norm_list[index] = alt_temp
            smoothing_list[index] = lambda_val
    if plot_S_curve:  
        figure('using NNLS output norm')
        rnorm_axis = array(rnorm_list)
        smoothing_axis = array(smoothing_list)
        plot(log10(smoothing_axis**2),rnorm_axis)
        show()
        figure();title('using LV norm')
        altnorm_axis = array(alt_norm_list)
        smoothing_axis = array(smoothing_list)
        plot(log10(smoothing_axis**2),altnorm_axis,'-.',c='k')
        xlabel(r'log($\alpha$)')
        ylabel(r'$\chi$($\alpha$)')
        gridandtick(gca())
        show() 
    if S_curve_guess:
        heel = raw_input("heel of S curve: ")
        heel_alpha = 10**heel
        heel_lambda = sqrt(heel_alpha)
        print "Alpha",heel_alpha
        print "Lambda",heel_lambda
        guess_lambda = heel_lambda


# In[ ]:


guess_lambda = 0.1 # Set to any desired variable
alpha = guess_lambda**2
print("Alpha",alpha)
print("Lambda",guess_lambda)


# In[ ]:


print("Estimating solution for guessed smoothing parameter...")
opt_vec = nnls_reg(guess_lambda)
solution = reshape(opt_vec,(Nx,Ny))
figure()
title(r'Est F(log$(p_{\frac{1}{2}}$),log$(T_{2})$, $\lambda$ = %0.3f'%(guess_lambda))
image(solution);show()


# In[ ]:


opt_val = mod_BRD(guess=guess_lambda,maxiter=20)
print("OPTIMIZED LAMBDA:",opt_val)


# In[ ]:


opt_vec = nnls_reg(opt_val)
solution = reshape(opt_vec,(Nx,Ny))


# In[ ]:


figure()
title(r'Est F(p$_{.5}$,log$(T_{2})$, $\lambda$ = %0.2f'%(opt_val))
image(solution);show()
    


# In[ ]:


nd_solution = nddata(solution,['pp5','log(T2)'])


# In[ ]:


nd_solution.setaxis('pp5',log_Nx_ax.copy())
nd_solution.setaxis('log(T2)',log_Ny_ax.copy())


# In[ ]:


figure();title(r' Estimated F(p$_{0.5}$,log($T_{2}$)), $\lambda$ = %0.6f'%opt_val)
nd_solution.contour(labels=False)
gcf().subplots_adjust(bottom=0.15)


# In[ ]:


10**-0.7


# In[ ]:




