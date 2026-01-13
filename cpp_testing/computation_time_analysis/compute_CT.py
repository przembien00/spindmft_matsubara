import argparse
import numpy as np

def InputProgramm():
    parser = argparse.ArgumentParser(description='All template options.')

    parser.add_argument('--samples', default=100000, type=int,
                    help='number of samples per core')

    parser.add_argument('--steps', default=500, type=int,
                    help='number of time steps')

    parser.add_argument('--iterations', default=4, type=int,
                    help='number of iterations')

    parser.add_argument('--stype', default="B", type=str,
                    help='symmetry type')
    return parser.parse_args()

args = InputProgramm()

m, dm, ps, dps = np.genfromtxt("sample_data.txt",unpack='True')
a1, da1, pt1, dpt1 = np.genfromtxt("steps_data_hs.txt",unpack='True')
a2, da2, pt2, dpt2 = np.genfromtxt("steps_data_l31.txt",unpack='True')

def computation_time(ps,a,pt,args):
    stype=0
    if args.stype == "C":
        stype=1
    elif args.stype == "D":
        stype=2
    time = (args.samples/1000)**ps * a[stype]*(args.steps)**pt[stype]
    return time*args.iterations

print("=====================================================================")
print("computation time for: M=", args.samples, ", steps=", args.steps, ",iterations=", args.iterations, ",stype=", args.stype)
print("at home : ", np.around(computation_time(ps[0],a1,pt1,args),2), "seconds")
print("at l31  : ", np.around(computation_time(ps[1],a2,pt2,args),2), "seconds")
print("=====================================================================")
