"""
Model class for simulation
"""
import numpy as np
import random as rd
import inspect
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from tqdm import tqdm
from scipy.sparse import csr_matrix
from multiprocessing import cpu_count
from multiprocessing import Pool
import os


class MG():
    """
    Model class intended to simulate
    """
    def __init__(self,
                 n_sim,
                 T=20,
                 BR=1,
                 BL=0,
                 I_1=1/3,
                 I_2=2/3,
                 delta_t=10**-3,
                 theta=1,
                 N=20,
                 fact=10,
                 start='full',
                 D=1,
                 mu = 10,
                 lambda_0 = 10000):
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        values.pop('self')

        for arg, val in values.items():
            setattr(self, arg, val)

        self.h = (self.BR-self.I_1)/self.N
        self.delta_x = self.h/self.fact
        self.k_C = round((self.BR-self.I_2)/self.h)
        self.k_I = round((self.I_2-self.I_1)/self.h)
        self.k = self.k_C+self.k_I
        self.k_pde = round((self.I_1-self.BL)/self.delta_x)
        self.k_pi = round((self.I_2-self.BL)/self.delta_x)

        # SSA Matrix
        Vr = np.eye(self.k, self.k)-np.eye(self.k, self.k, k=1)
        Vl = -np.eye(self.k, self.k-1)+np.eye(self.k, self.k-1, k=-1)

        self.V = np.concatenate((Vl, -Vr), axis=1)

        self.DC = np.array([self.D_C(self.I_1+i*self.h) for i in range(self.k)])/self.h**2
        self.muC = np.array([self.mu_C(self.I_1+(i+0.5)*self.h) for i in range(self.k)])

        # PDE Matrix
        self.Dpde = np.array([self.D-self.D_C(self.BL+i*self.delta_x) for i in range(self.k_pi+1)])/self.delta_x**2
        self.mu_pde = sp.diags(np.array([self.mu-self.mu_C(self.BL+(i+0.5)*self.delta_x) for i in range(self.k_pi)]))
        A = -sp.diags(self.Dpde[:-1]+self.Dpde[1:])+sp.diags(self.Dpde[1:-1], 1)+sp.diags(self.Dpde[1:-1], -1)
        A[0, 0] = -self.Dpde[0]
        A[-1, -1] = -self.Dpde[-2]
        self.A_1 = sp.eye(self.k_pi)-self.delta_t*self.theta*A
        self.A_2 = sp.eye(self.k_pi)+self.delta_t*(1-self.theta)*A

        self.B1 = sp.eye(self.k_pi)+self.delta_t*self.mu_pde-self.delta_t*A

        self.Source = np.array([self.lambda_0]+[0]*(self.k_pi-1))

    def D_C(self, x):
        if self.BL <= x <= self.I_1:
            return 0
        elif self.I_1 < x < self.I_2:
            return self.D*(x-self.I_1)/(self.I_2-self.I_1)
        elif self.I_2 <= x <= self.BR:
            return self.D

    def mu_C(self, x):
        if self.BL <= x <= self.I_1:
            return 0
        elif self.I_1 < x < self.I_2:
            return self.mu*(x-self.I_1)/(self.I_2-self.I_1)
        else:
            return self.mu

    def seuil(self, x):
        if 0 <= x < 1:
            return 0
        else:
            return x

    def simulate_one(self):
        # fonction qui retourne juste une matrice de simulation
        if self.start == 'empty':
            P = np.array([0]*(self.k_pi))
            S = np.array([0]*(self.k_I+self.k_C))
        elif self.start == 'full':
            P = np.array([33/self.h]*(self.k_pi))
            S = np.array([33]*(self.k_I+self.k_C))
        S = S.astype(float)
        P = P.astype(float)
        C = [np.concatenate((P, S/self.h))]
        for time_step in range(int(self.T/self.delta_t)):
            t = 0
            while t < self.delta_t:
                alphal = S*self.DC
                alphar = S[:-1]*self.DC[1:]
                beta = S*self.muC
                prop = np.concatenate((alphar, alphal, beta))
                a0 = sum(prop)
                if a0 == 0:
                    tau = self.delta_t+1
                else:
                    tau = -np.log(rd.random())/a0
                t = t+tau
                if t < self.delta_t:
                    j = np.random.multinomial(1, np.array(prop)/a0).nonzero()[0][0]
                    if j < self.k-1:              # jump to the right
                        x = np.array([1, S[j]]).min()
                        if j < self.k_I-1:         # jump in I
                            P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact] = P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact]-x/self.h
                            P[self.k_pde+(j+1)*self.fact:self.k_pde+(j+2)*self.fact] = P[self.k_pde+(j+1)*self.fact:self.k_pde+(j+2)*self.fact]+x/self.h
                        elif j == self.k_I-1:     # jump over I_2
                            P[self.k_pde+j*self.fact:] = P[self.k_pde+j*self.fact:]-x/self.h
                        S = np.add(S, x*self.V[:, j])
                    elif self.k-1 <= j < 2*self.k-1:                  # jump to the left
                        j = j-(self.k-1)
                        x = np.array([1,S[j]]).min()
                        if 0 < j < self.k_I:        # jump in I
                            P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact] = P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact]-x/self.h
                            P[self.k_pde+(j-1)*self.fact:self.k_pde+j*self.fact] = P[self.k_pde+(j-1)*self.fact:self.k_pde+j*self.fact]+x/self.h
                        elif j == self.k_I:        # jump over I_2
                            P[self.k_pde+(j-1)*self.fact:] = P[self.k_pde+(j-1)*self.fact:]+x/self.h
                        elif j == 0:          # jump left over I_1
                            P[self.k_pde:self.k_pde+self.fact] = P[self.k_pde:self.k_pde+self.fact]-x/self.h
                            P[self.k_pde-self.fact:self.k_pde] = P[self.k_pde-self.fact:self.k_pde]+x/self.h
                        S = np.add(S, x*self.V[:, j+self.k-1])
                    else:
                        j = j-(2*self.k-1)
                        x = np.array([1, S[j]]).min()
                        if j < self.k_I:
                            S[j] -= x
                            P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact] = P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact]-x/self.h
                        else:
                            S[j] -= x

            # P = spl.spsolve(A_1,A_2*P)-delta_t*mu_pde*P+D*delta_t*Source/delta_x
            # P = spl.spsolve(A_1+delta_t*mu_pde,P+delta_t*Source*D/h)
            # P = spl.spsolve(A_1,P-delta_t*mu_pde*P+delta_t*Source*D/delta_x)
            # P = A_2*P-delta_t*mu_pde*P+D*delta_t*Source/delta_x
            P = spl.spsolve(self.B1, P) + self.D*self.delta_t*self.Source/self.delta_x

            S[:self.k_I] = np.array([sum(P[self.k_pde+i*self.fact:self.k_pde+(i+1)*self.fact])*self.delta_x for i in range(self.k_I)])

            C.append(np.concatenate((P, S/self.h)))  # recording
        return np.array(C)

    def add_simulations(self,
                        C1,
                        C2):
        C = np.add(C1, C2)
        return C

    def simulate_n(self,
                   n):
        C = np.zeros((int(round(self.T/self.delta_t))+1, self.k_pi+self.k_I+self.k_C))  # mettre la matrice de zeros aux bonnes dimensions
        for i in tqdm(range(n)):
            C1 = self.simulate_one()
            C = self.add_simulations(C, C1)
        return C

    def multiprocess_sims(self,
                          N):
        nb_cpu = cpu_count()
        n_sims = N // nb_cpu
        args = [n_sims for i in range(nb_cpu)]
        pool = Pool(nb_cpu)
        results = pool.map(self.simulate_n, args)
        C = results[0]
        for i in tqdm(range(1, len(results))):
            C = self.add_simulations(C, results[i])
        self.mean = C/(n_sims*nb_cpu)
        return None

    def true_solution(self):
        n_pde = round((self.BR-self.BL)/self.delta_x)
        n1 = round(n_pde / 3)
        n2 = round(2 * n_pde / 3)
        M = np.zeros((n_pde, n_pde))
        M[-1][-1] = 1
        M[0][0] = 1
        A = (1+self.delta_t*self.mu)*np.eye(n_pde)-self.delta_t*(self.D/self.delta_x**2)*(np.eye(n_pde, k=1)+np.eye(n_pde, k=-1)-2*np.eye(n_pde)+M)
        S=np.array([self.delta_t*self.lambda_0*self.D/self.delta_x]+[0]*(n_pde-1))
        X2 = 0
        if self.start == 'empty':
            X = np.array([0]*n_pde)
            tab = [[0, 0, 0]]
        if self.start == 'full':
            X = np.array([990]*n_pde)
            tab = [[330, 330, 330]]
        for i in range(int(self.T/self.delta_t)):
            X = np.linalg.solve(A, X) + S
            tab.append([sum(X[:n1])*self.delta_x, sum(X[n1:n2])*self.delta_x, sum(X[n2:])*self.delta_x])
            if i == int(1/self.delta_t):
                X2 = X
        return tab, X, X2

    def save(self):
        if not os.path.isdir('experiment'):
            os.makedirs('experiment')
        if self.start == 'empty':
            np.save('experiment/mg_empty', self.mean)
        elif self.start == 'full':
            np.save('experiment/TEST_PROBLEM_3', self.mean)

    def load(self):
        if self.start == 'empty':
            self.mean = np.load('experiment/mg_empty.npy')
        elif self.start == 'full':
            self.mean = np.load('experiment/TEST_PROBLEM_3.npy')

    def plot(self):
        X_time = np.linspace(self.delta_t, self.T, int(self.T/self.delta_t))
        mean = self.mean
        tab, X, X2 = self.true_solution()
        if self.T > 1:
            time_plot = 1
            t_step = int(time_plot/self.delta_t)
            plt.figure(1)
            if self.start == 'empty':
                def u(x, ti):
                    return self.lambda_0*np.sqrt(self.D/self.mu)*np.cosh(np.sqrt(self.mu/self.D)*(x-1))/np.sinh(2*np.sqrt(self.mu/self.D))-self.lambda_0*self.D*np.exp(-self.mu*ti)/2/self.mu-4*self.lambda_0*self.D*sum(list(map(lambda y: np.cos(y*np.pi*(x+1)/2)*np.exp(-((self.D*(y*np.pi)**2/4)+self.mu)*ti)/(self.D*(y*np.pi)**2+4*self.mu), np.linspace(1, 1000, 1000))))
                # plt.plot(np.linspace(-1, 1, 1000), u(np.linspace(-1, 1, 1000), time_plot), 'k' + '--', linewidth=3)
            if self.start == 'full':
                plt.plot(np.linspace(self.BL, self.BR, len(X2)), X2,'k' + '--', linewidth=3)
                plt.bar(np.linspace(
                                    self.I_1,
                                    self.BR-self.h,
                                    int((self.BR-self.I_1)/self.h)
                                    ),
                        mean[t_step][self.k_pi:],
                        self.h,
                        align='edge')
                plt.plot(np.array([self.BL + (i+0.5)*self.delta_x for i in range(self.k_pi)]),
                         mean[t_step][:self.k_pi],
                         'r',
                         linewidth=3)
                # plt.title('{} simulations T={}'.format(self.n_sim, time_plot))
                plt.axis([0,1,0,4000])
                plt.title('T1')
                plt.plot()

        time_plot = self.T
        t_step = int(time_plot/self.delta_t)
        plt.figure(2)
        if self.start == 'empty':
            plt.plot(np.linspace(-1, 1, 1000),
                     u(np.linspace(-1, 1, 1000), time_plot),
                     'k' + '--',
                     linewidth=3)
        if self.start == 'full':
            plt.plot(np.linspace(self.BL, self.BR, len(X)), X,'k' + '--', linewidth=3)

        plt.bar(np.linspace(
                            self.I_1,
                            self.BR-self.h,
                            int((self.BR-self.I_1)/self.h)
                            ),
                mean[t_step][self.k_pi:],
                self.h,
                align='edge')
        plt.plot(np.array([self.BL+(i+0.5)*self.delta_x for i in range(self.k_pi)]),
                 mean[t_step][:self.k_pi],
                 'r',
                 linewidth=3)
        # plt.title('{} simulations T={}'.format(self.n_sim, time_plot))
        plt.axis([0,1,0,4000])
        plt.title('Tf')
        plt.plot()

        plt.figure(3)
        if self.start == 'empty':
            plt.plot(X_time, np.array([(sum(mean[t_time+1][self.k_pi+self.k_I:])*self.h-0)/1 for t_time in range(len(X_time))]),color='b', label='Compartment')
            plt.plot(X_time, np.array([(sum(mean[t_time+1][self.k_pde:self.k_pi])*self.delta_x-0)/1 for t_time in range(len(X_time))]),color='r', label='blending region')
            plt.plot(X_time, np.array([(sum(mean[t_time+1][:self.k_pde])*self.delta_x-0)/1 for t_time in range(len(X_time))]), color='g', label='pde')
        if self.start == 'full':
            plt.plot(X_time, np.array([(sum(mean[t_time+1][self.k_pi+self.k_I:])*self.h-tab[t_time+1][2])/tab[t_time+1][2] for t_time in range(len(X_time))]),color='b', label='Compartment')
            plt.plot(X_time, np.array([(sum(mean[t_time+1][self.k_pde:self.k_pi])*self.delta_x-tab[t_time+1][1])/tab[t_time+1][1] for t_time in range(len(X_time))]),color='r', label='blending region')
            plt.plot(X_time, np.array([(sum(mean[t_time+1][:self.k_pde])*self.delta_x-tab[t_time+1][0])/tab[t_time+1][0] for t_time in range(len(X_time))]), color='g', label='pde')
        plt.legend()
        plt.title('relative mass error')
        plt.plot()

        plt.show()
