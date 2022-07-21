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

class PureDiffusion():
    """
    Model class intended to simulate
    """
    def __init__(self,
                 n_sim,
                 L=1,
                 I_1=1/3,
                 I_2=2/3,
                 delta_t=10**-4,
                 theta=1,
                 T=1,
                 N=20,
                 delta_x=(2/3)/20/10,
                 N_particules=1000,
                 start='left',
                 D=1):
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        values.pop('self')

        for arg, val in values.items():
            setattr(self, arg, val)

        self.h = (self.L-self.I_1)/self.N
        self.fact = int(self.h/self.delta_x)
        self.k_C = int((self.L-self.I_2)/self.h)
        self.k_I = round((self.I_2-self.I_1)/self.h)
        self.k = self.k_C+self.k_I
        self.k_pde = int(self.I_1/self.delta_x)
        self.k_pi = int(self.I_2/self.delta_x)

        # SSA Matrix
        Vr = np.eye(self.k, self.k)-np.eye(self.k, self.k, k=1)
        Vl = -np.eye(self.k, self.k-1)+np.eye(self.k, self.k-1, k=-1)
        self.V = np.concatenate((Vl, -Vr), axis=1)
        self.DC = np.array([self.D_C(self.I_1+i*self.h) for i in range(self.k_C+self.k_I)])/self.h**2

        # PDE Matrix
        Dpde = np.array([self.D-self.D_C(i*self.delta_x) for i in range(self.k_pi+1)])/self.delta_x**2
        A = -sp.diags(Dpde[:-1]+Dpde[1:])+sp.diags(Dpde[1:-1], 1)+sp.diags(Dpde[1:-1], -1)
        A[0, 0] = -Dpde[0]
        A[-1, -1] = -Dpde[-2]
        self.A_1 = sp.eye(self.k_pi)-self.delta_t*self.theta*A
        self.A_2 = sp.eye(self.k_pi)+self.delta_t*(1-self.theta)*A

        # create self.MEAN: # TODO: mettre les self
        # self.MEAN = np.zeros((int(round(self.T/self.delta_t))+1,
        #                      self.k_pi+self.k_I+self.k_C))

    def D_C(self, x):
        if 0 <= x <= self.I_1:
            return 0
        elif self.I_1 < x < self.I_2:
            return self.D*(x-self.I_1)/(self.I_2-self.I_1)
        elif self.I_2 <= x <= self.L:
            return self.D

    def seuil(self, x):
        if 0 <= x < 1:
            return 0
        else:
            return x

    def simulate_one(self):
        # fonction qui retourne juste une matrice de simulation
        if self.start == 'left':
            P = np.concatenate((np.array([self.N_particules/self.I_1]*self.k_pde),
                                [0]*(self.k_pi-self.k_pde)))
            S = np.array([0]*self.N)
        # elif self.start == 'right':
        #     P = np.array([0]*(self.k_pi))
        #     S = np.concatenate(([0]*self.k_I, [float(self.N_particules/self.k_C)]*self.k_C))
        elif self.start == 'equilibrium':
            a = self.N_particules*(self.I_2-0)/self.L
            b = self.N_particules*(self.L-self.I_1)/self.L
            P = np.array([a/self.I_2]*self.k_pi)
            S = np.array([b/self.N]*self.N)
        S = S.astype(float)
        P = P.astype(float)
        C = [np.concatenate((P, S/self.h))]
        for time_step in range(int(self.T/self.delta_t)):
            t = 0
            while t < self.delta_t:
                alphal = S*self.DC
                alphar = S[:-1]*self.DC[1:]
                prop = np.concatenate((alphar, alphal))
                a0 = sum(prop)
                if a0 == 0:
                    tau = self.delta_t+1
                else:
                    tau = -np.log(rd.random())/a0
                t = t+tau
                if t <= self.delta_t:
                    j =\
                        np.random.multinomial(1, np.array(prop)/a0).nonzero()[0][0]
                    # which jump take place
                    if j < self.k-1:               # jump to the right
                        if j < self.k_I-1:         # jump in I
                            P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact] = P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact]-1/self.h
                            P[self.k_pde+(j+1)*self.fact:self.k_pde+(j+2)*self.fact] =\
                                P[self.k_pde+(j+1)*self.fact:self.k_pde+(j+2)*self.fact]+1/self.h
                        elif j == self.k_I-1:      # jump over I_2
                            P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact] =\
                                P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact]-1/self.h
                        S = np.add(S, self.V[:, j])
                    else:                          # jump to the left
                        j = j-(self.k-1)
                        if 0 < j < self.k_I:            # jump in I
                            P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact] =\
                                P[self.k_pde+j*self.fact:self.k_pde+(j+1)*self.fact]-1/self.h
                            P[self.k_pde+(j-1)*self.fact:self.k_pde+j*self.fact] =\
                                P[self.k_pde+(j-1)*self.fact:self.k_pde+j*self.fact]+1/self.h
                        elif j == self.k_I:            # jump over I_2
                            P[self.k_pde+(j-1)*self.fact:self.k_pde+j*self.fact] =\
                                P[self.k_pde+(j-1)*self.fact:self.k_pde+j*self.fact]+1/self.h
                        elif j == 0:              # jump left over I_1
                            P[self.k_pde:self.k_pde+self.fact] = P[self.k_pde:self.k_pde+self.fact]-1/self.h
                            P[self.k_pde-self.fact:self.k_pde] = P[self.k_pde-self.fact:self.k_pde]+1/self.h
                        S = np.add(S, self.V[:, j+ self.k-1])

            P = spl.spsolve(self.A_1, self.A_2*P)  # PDE update
            S[:self.k_I] =\
                np.array([self.seuil(
                                self.delta_x * sum(P[self.k_pde+i*self.fact:self.k_pde+(i+1) * self.fact])
                                ) for i in range(self.k_I)])
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
        args[0] += N % nb_cpu
        pool = Pool(nb_cpu)
        results = pool.map(self.simulate_n, args)
        C = results[0]
        for i in tqdm(range(1, len(results))):
            C = self.add_simulations(C, results[i])
        self.mean = C / N
        return None

    def left_true_solution(self):
        n_pde = int(self.L/self.delta_x)
        n1 = int(n_pde /3)
        n2 = int(2 * n_pde /3)
        M = np.zeros((n_pde, n_pde))
        M[-1][-1] = 1
        M[0][0] = 1
        A = np.eye(n_pde)-self.delta_t*(self.D/self.delta_x**2)*(np.eye(n_pde, k=1)+np.eye(n_pde, k=-1)-2*np.eye(n_pde)+M)
        X = np.array([3]*n1+[0]*(n_pde - n1))
        tab = [[1, 0, 0]]
        X2 = 0
        for i in range(int(self.T/self.delta_t)):
            X = np.linalg.solve(A, X)
            tab.append([sum(X[:n1])*self.delta_x, sum(X[n1:n2])*self.delta_x, sum(X[n2:])*self.delta_x])
            if i == int(0.1/self.delta_t):
                X2 = X
        return tab, X, X2

    def save(self):
        if not os.path.isdir('experiment'):
            os.makedirs('experiment')
        if self.start == 'left':
            np.save('experiment/TEST_PROBLEM_2', self.mean)
        elif self.start == 'equilibrium':
            np.save('experiment/TEST_PROBLEM_1', self.mean)

    def load(self):
        if self.start == 'left':
            self.mean = np.load('experiment/TEST_PROBLEM_2.npy')
        elif self.start == 'equilibrium':
            self.mean = np.load('experiment/TEST_PROBLEM_1.npy')

    def plot(self):
        X_time = np.linspace(self.delta_t, self.T, int(self.T/self.delta_t))
        mean = self.mean
        time_plot = self.T
        t_step = int(time_plot/self.delta_t)

        plt.figure(1)
        if self.start == 'left':
            tab, x1, x2 = self.left_true_solution()
            plt.plot(np.linspace(0, 1, len(x1)),
                     x1,
                     'k' + '--',
                     linewidth=3)
        if self.start == 'equilibrium':
            plt.plot(np.linspace(0, 1, 1000),
                     np.ones(1000),
                     'k' + '--',
                     linewidth=3)
        plt.bar(np.linspace(
                            self.I_1,
                            self.L-self.h,
                            self.N),
                mean[t_step][self.k_pi:] * 0.001,
                self.h,
                align='edge')
        plt.plot(np.array([(i+0.5)*self.delta_x for i in range(self.k_pi)]),
                 mean[t_step][:self.k_pi]*0.001,
                 'r',
                 linewidth=3)
        # plt.title('{} simulations T={}'.format(self.n_sim, time_plot))
        plt.axis([0,1,0,3.5])
        plt.title('T={}'.format(time_plot))
        plt.plot()

        if self.T > 0.15:
            time_plot = 0.1
            t_step = int(time_plot/self.delta_t)
            plt.figure(2)
            if self.start == 'left':
                plt.plot(np.linspace(0, 1, len(x1)),
                         x2,
                         'k' + '--',
                         linewidth=3)
            if self.start == 'equilibrium':
                plt.plot(np.linspace(0, 1, 1000),
                         np.ones(1000),
                         'k' + '--',
                         linewidth=3)
            plt.bar(np.linspace(
                                self.I_1,
                                self.L-self.h,
                                self.N),
                    mean[t_step][self.k_pi:] * 0.001,
                    self.h,
                    align='edge')
            plt.plot(np.array([(i+0.5)*self.delta_x for i in range(self.k_pi)]),
                     mean[t_step][:self.k_pi]*0.001,
                     'r',
                     linewidth=3)
            plt.axis([0,1,0,3.5])
            plt.title('T = 0.1')
            plt.plot()

        plt.figure(3)
        if self.start == 'equilibrium':
            n_equi_split = self.N_particules / 3
            plt.plot(X_time, np.array([(sum(mean[t_time+1][self.k_pi+self.k_I:])*self.h-n_equi_split)/n_equi_split for t_time in range(len(X_time))]),color='b', label='Compartment')
            plt.plot(X_time, np.array([(sum(mean[t_time+1][self.k_pde:self.k_pi])*self.delta_x-n_equi_split)/n_equi_split for t_time in range(len(X_time))]),color='r', label='blending region')
            plt.plot(X_time, np.array([(sum(mean[t_time+1][:self.k_pde])*self.delta_x-n_equi_split)/n_equi_split for t_time in range(len(X_time))]), color='g', label='pde')
        if self.start == 'left':
            plt.plot(X_time, np.array([(sum(mean[t_time+1][self.k_pi+self.k_I:])*self.h-tab[t_time+1][2]*1000)/1000/tab[t_time+1][2] for t_time in range(len(X_time))]),color='b', label='Compartment')
            plt.plot(X_time, np.array([(sum(mean[t_time+1][self.k_pde:self.k_pi])*self.delta_x-tab[t_time+1][1]*1000)/1000/tab[t_time+1][1] for t_time in range(len(X_time))]),color='r', label='blending region')
            plt.plot(X_time, np.array([(sum(mean[t_time+1][:self.k_pde])*self.delta_x-tab[t_time+1][0]*1000)/1000/tab[t_time+1][0] for t_time in range(len(X_time))]), color='g', label='pde')
        plt.axis([0, self.T,-0.02, 0.02])
        plt.legend()
        plt.title('relative mass error')
        plt.plot()

        # if self.start == 'left':
        #     plt.figure(4)
        #     plt.plot(X_time, np.array([(sum(mean[t_time+1][self.k_pi+self.k_I:])*self.h-tab[t_time+1][2]*1000)/1000 for t_time in range(len(X_time))]),color='b', label='Compartment')
        #     plt.plot(X_time, np.array([(sum(mean[t_time+1][self.k_pde:self.k_pi])*self.delta_x-tab[t_time+1][1]*1000)/1000 for t_time in range(len(X_time))]),color='r', label='blending region')
        #     plt.plot(X_time, np.array([(sum(mean[t_time+1][:self.k_pde])*self.delta_x-tab[t_time+1][0]*1000)/1000 for t_time in range(len(X_time))]), color='g', label='pde')
        #     plt.legend()
        #     plt.title('mass error normalized')
        #     plt.plot()
        plt.show()
