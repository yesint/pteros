from pteros import *
import numpy as np
from sklearn.cluster import AgglomerativeClustering,Birch

#----------------------------------------------------
def do_clustering(m,Ncl,output_sel):
    # Agglomerative clustering
    clustering = AgglomerativeClustering(affinity="precomputed",n_clusters=Ncl,linkage='complete').fit(m)

    print(clustering.labels_)

    # For each cluster find the member with least average difference
    cl = {l:[] for l in range(Ncl)}

    for i,l in enumerate(clustering.labels_):
        cl[l].append(i)

    for l in range(Ncl):    
        print(f'cluster {l}: {cl[l]}')
        min_el=-1
        min_val=1e6
        for el1 in cl[l]:
            ave = 0.0
            for el2 in cl[l]:
                ave += m[el2,el1]
            ave /= len(cl[l])
            
            if ave<min_val:
                min_val = ave
                min_el = el1
                
            print(f' el {el1} -> {ave}')
            
        print(f'Central element: {min_el}')
        
        output_sel.write(f'cluster_{l}.pdb',min_el,min_el)

#----------------------------------------------------


print('Loading trajectory...')
s=System()
s.set_filter('not resname SOL NA CL')
s.load('topol.tpr')
s.frame_delete(0,0) # We don't need first frame from topol.tpr
s.load('traj_comp.xtc',2500,-1)

fit_sel = s('resid 75 76 77 78 80 83 85 87 141 142 143 177 191 193 206 365 366 367 368 369 370 391 392')
out_sel = s()

print(f'RMSD selection size: {fit_sel.size()}')

print('Unwrapping selection...')
for fr in range(s.num_frames()):
    out_sel.set_frame(fr)
    out_sel.unwrap_bonds(d=0.25)

print('Computing RMSD matrix...')
m = rmsd_matrix(fit_sel)

print('Performing clustering...')
do_clustering(m,200,out_sel)



