import numpy as np 
import matplotlib.pyplot as plt 
from itertools import combinations
import vectors as v
import copy

'''
todo's:
incorporate spin and lattice classes everywhere
write class for distance vectors?
'''

# rounding routines:
r2 = lambda value : np.around(value,2) # round values for terminal output
r3 = lambda value : np.around(value,3) # round values for terminal output
r4 = lambda value : np.around(value,4) # round values for terminal output
r5 = lambda value : np.around(value,5) # round values for terminal output

origin = np.array([0,0,0])

# =========================================
# ================ Classes ================
# =========================================
class Spin:
    r = origin
    basis_index = 0
    Bravais_index = [0,0,0]
    sort_index = 0

    def __init__(self,r,basis_index,Bravais_index,sort_index=0):
        self.r = r 
        self.basis_index = basis_index
        self.Bravais_index = Bravais_index
        self.sort_index = sort_index

    def __eq__(self,other):
        if isinstance(other,Spin):
            return (self.basis_index == other.basis_index and
                    self.Bravais_index == other.Bravais_index)
        return False


class Spin_Lattice_3D:
    lattice = []
    lattice_vectors = None
    basis_vectors = None
    sort_indices = None
    central_spin_index = None

    def __init__(self,num_cells_in_dir,lattice_vectors,basis_vectors=[origin],sort_indices="single",exclude=[]):
        if sort_indices == "single": # just a single spin sort
            sort_indices = [0 for i in range(len(basis_vectors))]
        linear_index = 0
        for h in range(-num_cells_in_dir,num_cells_in_dir+1):
            for k in range(-num_cells_in_dir,num_cells_in_dir+1):
                for l in range(-num_cells_in_dir,num_cells_in_dir+1):
                    for basis_index, (sort_index,b) in enumerate(zip(sort_indices,basis_vectors)):
                        if not np.any( [excl_bi == basis_index and excl_Bis == [h,k,l] for excl_bi, excl_Bis in exclude] ): # exclude some spins
                            pos = h*lattice_vectors[0] + k*lattice_vectors[1] + l*lattice_vectors[2] + b        
                            self.lattice.append(Spin(pos,basis_index,[h,k,l],sort_index=sort_index))
                            if basis_index == 0 and [h,k,l] == [0,0,0]:
                                self.central_spin_index = linear_index
                            linear_index += 1
        self.lattice_vectors = lattice_vectors
        self.basis_vectors = basis_vectors
        self.sort_indices = sort_indices
    
    def find_from_index(self,basis_index,Bravais_index):
        for spin in self.lattice:
            if spin.basis_index == basis_index and spin.Bravais_index == Bravais_index:
                return True 
        return False

    def find_from_position(self,r):
        for spin in self.lattice:
            if np.array_equal(spin.r,r):
                return True 
        return False

    def return_from_index(self,basis_index,Bravais_index):
        for spin in self.lattice:
            if spin.basis_index == basis_index and spin.Bravais_index == Bravais_index:
                return spin 
        raise RuntimeError( "spin does not exist in lattice" )

    def return_from_position(self,r):
        for spin in self.lattice:
            if np.array_equal(spin.r,r):
                return spin 
        raise RuntimeError( "spin does not exist in lattice" )

    def return_central_spin(self):
        if self.central_spin_index != None:
            return self.lattice[self.central_spin_index]
        raise RuntimeError( "spin does not exist in lattice" )
        

# =========================================
# ================ Lattices ===============
# =========================================
# creates a lattice from lattice vectors (spanning the unit cell) and basis vectors
# the lattice contains 3*(2*num_cells_in_dir+1) cells
# some spins can be excluded manually
def create_3D_lattice(lattice_vectors,basis_vectors,num_cells_in_dir,excluded_spins=[]):
    # create lattice:
    lattice = [] 
    spin_types = []
    for h in range(-num_cells_in_dir,num_cells_in_dir+1):
        for k in range(-num_cells_in_dir,num_cells_in_dir+1):
            for l in range(-num_cells_in_dir,num_cells_in_dir+1):
                for b_index, b in enumerate(basis_vectors):
                    r = h*lattice_vectors[0] + k*lattice_vectors[1] + l*lattice_vectors[2] + b
                    lattice.append(r)
                    spin_types.append(b_index)
    return exclude2(lattice,spin_types,excluded_spins)

# exclude spins from a lattice
def exclude(lattice_,excluded_spins):
    lattice = copy.deepcopy(lattice_)
    for excluded_spin in excluded_spins:
        lattice = [spin for spin in lattice if not np.array_equal(spin, excluded_spin)]
    return lattice

def exclude2(lattice_,spin_types_,excluded_spins):
    lattice, spin_types = copy.deepcopy(lattice_), copy.deepcopy(spin_types_)
    for excluded_spin in excluded_spins:
        (lattice, spin_types) = np.transpose([(spin, spin_type) for spin, spin_type in zip(lattice, spin_types) if not np.array_equal(spin, excluded_spin)])
    return lattice, spin_types

# creates a lattice from lattice vectors (spanning the unit cell) and basis vectors
# the lattice is build around the rhomboid defined by the given Bravais numbers [hkl] of two cells
# the thickness of the coeat is handed over as parameter
# some spins can be excluded manually
def create_3D_coat(lattice_vectors,basis_vectors,coat_thickness,cell_i,cell_j,excluded_spins):
    # create ranges:
    ranges = []
    for dir in range(3):
        min = np.min([cell_i[dir],cell_j[dir]]) 
        max = np.max([cell_i[dir],cell_j[dir]])
        ranges.append( range(min-coat_thickness,max+coat_thickness+1) )
    # create lattice:
    lattice = [] 
    spin_types = []
    for h in ranges[0]:
        for k in ranges[1]:
            for l in ranges[2]:
                for b_index, b in enumerate(basis_vectors):
                    r = h*lattice_vectors[0] + k*lattice_vectors[1] + l*lattice_vectors[2] + b
                    lattice.append(r)
                    spin_types.append(b_index)
    return exclude2(lattice,spin_types,excluded_spins)

# determine the Bravais lattice indices [hkl] for a given site
def determine_Bravais_numbers(ri,lattice_vectors,basis_vectors,max_number=5):
    for h in range(-max_number,max_number+1):
        for k in range(-max_number,max_number+1):
            for l in range(-max_number,max_number+1):
                for b_index, b in enumerate(basis_vectors):
                    rj = h*lattice_vectors[0] + k*lattice_vectors[1] + l*lattice_vectors[2] + b
                    if np.array_equal(ri,rj):
                        return [h,k,l], b_index


# =================================================
# ============== Shapes and Routes ================
# =================================================
# finds all possible combinations of placing k entities in N spots (numbered from 0 to N-1)
# the entities are assumed indistinguishable so that any permutations are thrown out
def generate_sorted_lists(k, N):
    elements = list(range(N))
    all_combinations = combinations(elements, k)
    result = [list(combination) for combination in all_combinations]
    return result

# creates all combinations of shapes assuming a couple of nodes fixed (gen_shape)
# and num_remaining variable nodes which can be placed anywhere in free_spots
def create_all_combinations_of_shapes(gen_shape,free_spots,num_remaining):
    shapes = []
    possible_index_ranges = generate_sorted_lists(num_remaining, len(free_spots))
    for possible_index_range in possible_index_ranges:
        new_shape = gen_shape.copy()
        for index in possible_index_range:
            new_shape.append(free_spots[index])
        shapes.append(new_shape)
    return shapes

# create all routes from ri and rj with n nodes, i.e., n+1 links
# any routes that are not within the allowed_spots are excluded
def create_all_routes(n,ri,type_i,rj,type_j,dvs_per_type,spts_per_type,allowed_spots="All"):
    if allowed_spots == "All":
        condition = lambda *args : True
    else:
        condition = lambda new_r : np.any( [np.array_equal(new_r,rk) for rk in allowed_spots] )
    routes = [[ri]] # list of lists of arrays (ri, d to next spin, d's to next spin, ...)
    current_types = [type_i] # list of current spin types corresponding to last entry of routes
    current_rs = [ri]
    for step in range(0,n+1): # create all possible routes from spin_i with n+1 links
        new_routes = []
        new_types = []
        new_rs = []
        for route_index, route in enumerate(routes):
            previous_type = current_types[route_index]
            previous_r = current_rs[route_index]
            for new_d, new_type in zip(dvs_per_type[previous_type],spts_per_type[previous_type]): # loop over all possibilities for the next link
                new_r = previous_r + new_d
                if condition( new_r ): # is the position of the new site allowed?
                    new_route = route + [new_d]
                    new_routes.append(new_route)
                    new_types.append(new_type)
                    new_rs.append(new_r)
        routes = new_routes
        current_types = new_types
        current_rs = new_rs
    valid_routes = []
    for route, type in zip(routes,current_types): # find all valid routes, i.e., those which end at spin_j
        tot_dv = np.sum(route, axis=0)
        if np.array_equal(tot_dv, rj):
            if type != type_j: # check whether the types match
                raise RuntimeError("Distance vector matches but the spin types not. Something is not right...")
            valid_routes.append(route[1:-1]) # leave out spin_i and spin_j => the route consists only of the distance vectors for the n+1 links
    return valid_routes


# add the spins on a route to a shape
def add_route_to_shape(shape,ri,links):
    rk = ri
    for d in links:
        rk = rk + d
        if not np.any( [np.array_equal(rk,rl) for rl in shape] ): # if spin not in shape
            shape.append(rk)
    return shape


# =================================================
# ================ Dipolar Couplings ==============
# =================================================
def I_of_theta(theta): # angular dependence of the couplings
    return (1-3*np.cos(theta)**2)

def I_of_costheta(costheta): # angular dependence of the couplings
    return (1-3*costheta**2)

def J_of_costheta(r,costheta): # couplings
    return (1-3*costheta**2)/r**3

def J_of_r_plane(r_vec):
    r = v.norm(r_vec)
    cos2phi = np.cos( 2*np.arccos(r_vec[0]/r) )
    return cos2phi/r**3 

def J_of_r_plane_tilted(r_vec,phi0):
    r = v.norm(r_vec)
    if r_vec[1]>0:
        phi = np.arccos(r_vec[0]/r)
    else:
        phi = -np.arccos(r_vec[0]/r) # +2pi can be ignored
    return np.cos( 2*(phi-phi0) )/r**3 

def J_of_r(r_vec,nB):
    r = v.norm(r_vec)
    costheta = v.scalar(r_vec,nB) / r
    return J_of_costheta(r, costheta)

def J_of_r_simple(r_vec):
    r = v.norm(r_vec)
    return 1/r**3

# =================================================
# ================ Coupling Constants =============
# =================================================
# scheme: J#1_p#2_of_#3 || #1 = order/name, #2 = power, #3 = dependence

def JQ_p2_of_r(positions,ri,*args,coupling_of_r=J_of_r):
    c = 0
    for rj in positions:
        if not np.array_equal(rj,ri):
            c += coupling_of_r(rj-ri,*args)**2
    return c

def JQ_p2_of_i(positions,i,*args,coupling_of_r=J_of_r):
    c = 0
    ri = positions[i]
    for j, rj in enumerate(positions):
        if j!=i:
            c += coupling_of_r(rj-ri,*args)**2
    return c

def JT_p4_of_r(positions,ri,*args,coupling_of_r=J_of_r):
    c = 0
    for rj in positions:
        if not np.array_equal(rj,ri):
            c += coupling_of_r(rj-ri,*args)**4
    return c

def JT_p4_of_i(positions,i,*args,coupling_of_r=J_of_r):
    c = 0
    ri = positions[i]
    for j, rj in enumerate(positions):
        if j!=i:
            c += coupling_of_r(rj-ri,*args)**4
    return c

# compute coupling constant for given links (list of relative distance vectors)
# the vectors are added upon ri 
# neither ri nor rj is included
def Jlinks_pN_of_r(ri,rj,links,*args,coupling_of_r=J_of_r):
    c = 1
    r = ri
    for r_rel in links:
        r_new = r + r_rel
        c *= coupling_of_r(r_new-r,*args)**2
        r = r_new
    return c * coupling_of_r(rj-r,*args)**2

# compute coupling constant for given nodes (list of positions)
# the nodes include neither ri nor rj
def Jnodes_pN_of_r(ri,rj,nodes,*args,coupling_of_r=J_of_r):
    c = 1
    r = ri
    for r_new in nodes:
        c *= coupling_of_r(r_new-r,*args)**2
        r = r_new
    return c * coupling_of_r(rj-r,*args)**2

def JC_p4_of_r(positions,ri,rj,*args,coupling_of_r=J_of_r):
    c = 0
    for rk in positions:
        if not np.array_equal(rk,ri) and not np.array_equal(rk,rj):
            c += Jnodes_pN_of_r(ri,rj,[rk],*args,coupling_of_r=coupling_of_r)
    return c

def JC_p4_of_i( positions, i, j,*args,coupling_of_r=J_of_r):
    return JC_p4_of_r(positions,positions[i],positions[j],*args,coupling_of_r=coupling_of_r)

# check whether a position is in a list of positions
def isin(r0,positions):
    for r in positions:
        if np.array_equal(r0,r):
            return True
    return False


# =================================================
# ================ Coordination Numbers ===========
# =================================================
# compute the effective coordination number for given positions (list of positions)
def zeff_of_r( positions, ri, *args,coupling_of_r=J_of_r):
    JQ_p4 = ( JQ_p2_of_r(positions,ri,*args,coupling_of_r=coupling_of_r) )**2
    JT_p4 = JT_p4_of_r(positions,ri,*args,coupling_of_r=coupling_of_r)
    return JQ_p4 / JT_p4

def zeff_of_i( positions, i, *args,coupling_of_r=J_of_r):
    JQ_p4 = ( JQ_p2_of_i(positions,i,*args,coupling_of_r=coupling_of_r) )**2
    JT_p4 = JT_p4_of_i(positions,i,*args,coupling_of_r=coupling_of_r)
    return JQ_p4 / JT_p4


# =================================================
# ===================== Others ====================
# =================================================
# sort the pair correlations (shells) according to the distances and couplings
def sort_shells(lattice,num_shells,*args,r0=origin,coupling_of_r=J_of_r):
    # compute couplings:
    distances, abs_Js, distancevectors, occurrences = [], [], [], []
    for r in lattice:
        dist = r5(v.norm(r))
        abs_J = r5(np.abs(coupling_of_r(r-r0,*args)))
        pos = np.where((distances == dist) & (abs_Js == abs_J))[0]
        if pos.size == 0: # add if pair (dist,J) does not exist yet
            distances.append(dist)
            abs_Js.append(abs_J)
            distancevectors.append(r)
            occurrences.append(1)
        else:
            occurrences[pos[0]] = occurrences[pos[0]]+1

    # first sort them in increasing order with the distance, then in descending order with the couplings:
    argsort = np.lexsort((-np.array(abs_Js),distances))
    distances_sorted, abs_Js_sorted, distancevectors_sorted, occurrences_sorted = [0.],[0.],[r0],[1] # add the origin automatically
    for sort, shell_index in zip(argsort,range(num_shells)):
        distances_sorted.append(distances[sort])
        abs_Js_sorted.append(abs_Js[sort])
        distancevectors_sorted.append(distancevectors[sort])
        occurrences_sorted.append(occurrences[sort])

    # return
    return distances_sorted, abs_Js_sorted, distancevectors_sorted, occurrences_sorted



# find the L most important couplings to a spin (in case of equal couplings L is automatically increased)
# the spin should be in the central unit cell of the lattice
def find_most_important_couplings(L,lattice_,spin_types_,ri,*args,coupling_of_r=J_of_r):
    # exclude spins:
    lattice, spin_types = exclude2(lattice_,spin_types_,[ri])

    # consider all couplings to ri:
    abs_Js, distancevectors, dst_types = [], [], []
    for rj, tj in zip(lattice,spin_types):
        abs_Js.append( r5(np.abs(coupling_of_r(rj-ri,*args))) )
        distancevectors.append(rj-ri)
        dst_types.append(tj)

    # sort them in descending order with the absolute couplings:
    argsort = np.flip(np.argsort(np.array(abs_Js)))
    distancevectors_sorted, abs_Js_sorted, spin_types_sorted = [], [], []
    n = 0
    for sort in argsort:
        abs_J = abs_Js[sort]
        if n<L or abs_J == abs_Js_sorted[-1]: # second condition ensures that no equally important couplings are thrown away
            abs_Js_sorted.append(abs_J)
            distancevectors_sorted.append(distancevectors[sort])
            spin_types_sorted.append(dst_types[sort])
        else:
            break
        n += 1

    return distancevectors_sorted, abs_Js_sorted, spin_types_sorted

# find the L most important couplings to all spins in the central unit cell
def find_most_important_couplings_for_all_types(L,lattice,spin_types,basis_vectors,*args,coupling_of_r=J_of_r):
    dvs_per_type, spts_per_type, aJs_per_type = [], [], []
    for type, bv in enumerate(basis_vectors):
        dvs, absJs, spts = find_most_important_couplings(L,lattice,spin_types,bv,*args,coupling_of_r=coupling_of_r)
        dvs_per_type.append(dvs)
        spts_per_type.append(spts)
        aJs_per_type.append(absJs)
    return dvs_per_type, spts_per_type, aJs_per_type
