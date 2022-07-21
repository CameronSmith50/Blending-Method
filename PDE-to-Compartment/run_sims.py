from pure_diffusion import PureDiffusion
from morphogen_gradient import MG

chemical_reaction = 'TEST_PROBLEM_3'   # 'TEST_PROBLEM_X'

if chemical_reaction == 'TEST_PROBLEM_1':
    n_sim = 500
    purediff = PureDiffusion(n_sim)
    purediff.multiprocess_sims(n_sim)
    purediff.save()

elif chemical_reaction == 'TEST_PROBLEM_2':
    n_sim = 500
    purediff = PureDiffusion(n_sim, start='equilibrium')
    purediff.multiprocess_sims(n_sim)
    purediff.save()

elif chemical_reaction == 'TEST_PROBLEM_3':
    n_sim = 100
    mg = MG(n_sim,T=4)
    mg.multiprocess_sims(n_sim)
    mg.save()
