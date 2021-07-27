import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.constants import hbar,k
import math
from random import randint
from numpy import random

#specify constants
l_size=int(30)
eps_1=1.41
eps_2=1
w=10**13
T=100
eps_1_2=-0.5
eps_1_1=-0.5
eps_2_2=-1
beta = 1/(k*T)	

#create and randomly fill 100 empty lattice with molecules in configuration 1 or 2
def create_lattice(lattice_size):
	x=np.zeros(lattice_size)
	y=np.zeros(lattice_size)
	x,y=np.meshgrid(x,y)
	chosenx=[]
	for _ in range(101):
		cx = randint(0, lattice_size-1)
		chosenx.append(cx)
	choseny=[]
	for _ in range(101):
		cy = randint(0, lattice_size-1)
		choseny.append(cy)	
	randlist=[]
	for _ in range(101):	
		n = random.random()
		randlist.append(n)
	for x1,y1,r in zip(chosenx,choseny,randlist):
		r=random.rand()
		if 0<r<0.7:
			x[x1,y1]=1	
		if r>=0.7:
			x[x1,y1]=2	

	return x		

#functiion used to create a lattice of vibrational energies where each entry corresponds to the vibrational energy of the molecule at the spcified point
def create_vib_lattice(l_size):
	E_x=np.zeros(l_size)
	E_y=np.zeros(l_size)
	E_x,E_y=np.meshgrid(E_x,E_y)
	
	return E_x
	

#Apply function to create vibrational lattice.
x=create_lattice(l_size)
E_vibr=create_vib_lattice(l_size)

plt.imshow(x)
plt.legend(x)
plt.show()

#assign a random virbational energy to every site occupied by a moleule in configuration 2.
def fill_vibr_lattice(E_vibr):
	for i in range(0,l_size):
			for k in range(0,l_size):
				if x[i][k]==2:
					E_vibr[i][k]=random.rand()*hbar*w
	return E_vibr				

E_vibr=fill_vibr_lattice(E_vibr)
	
#Pich a random point on the lattice(cpx,cpy) and check how many free neighbouts it has(j).
def has_free(cpx,cpy):
	R=[1,-1,0,0]
	L=[0,0,1,-1]
	x_n=[]
	y_n=[]
	j=0
	for a,b in zip(R,L):
		nx=cpx+a
		ny=cpy+b
		#periodic boundary condiitons
		if nx ==l_size:
			nx=0
		if nx==-1:
			nx=l_size-1
		if ny==l_size:
			ny=0
		if ny==-1:
			ny=l_size-1	
		#add to point to nearest neighbour list	
		if x[nx][ny]==0:
			j=j+1	
			
	
	
		
	return j
	
#function to find the x and y coordinates of free nearest neighbour site n, for a chosen point on the lattice (cpx),cpy) 	
def nearest_neighbours(cpx,cpy):
	R=[1,-1,0,0]
	L=[0,0,1,-1]
	x_n=[]
	y_n=[]
	for a,b in zip(R,L):
			nx=cpx+a
			ny=cpy+b
			#periodic boundary conditions
			if nx ==l_size:
				nx=0
			if nx==-1:
				nx=l_size-1
			if ny==l_size:
				ny=0
			if ny==-1:
				ny=l_size-1	
			#add site to list if it is free	
			if x[nx][ny]==0:	
				x_n.append(nx)
				y_n.append(ny)	
	
	return x_n, y_n			



#function used to calculate the contribution total binding energy. Count the number of site 1 and the number of site 2 and multiply by epsilon 1 and epsilon 2.
def calc_energy_site(x):
	n_1=0
	n_2=0
	total_E1=0
	total_E2=0
	for i in range(0,l_size):
		for k in range(0,l_size):
			if x[i][k]==1:
				n_1=n_1+1
			total_E1=n_1*eps_1
			if x[i][k]==2.0:
				n_2=n_2+1
			total_E2=n_2*eps_2
	#return total binding energy		
	return total_E1+total_E2
	
# function to calculate the total contributon of interaction energies to the total energies(I feel I maybe counting each interaction twice here)
def calc_enery_inter(x):
	R=[1,-1,0,0]
	L=[0,0,1,-1]
	for i in range(0,l_size):
		for k in range(0,l_size):
			#define the pairs of interaction, so n_1_1 is the interaction between 2 neighbouring molecules in configuration 1.	
			n_1_1=0
			n_1_2=0
			n_2_2=0
			n_1_0=0
			n_2_0=0
			#find the number of each type of interaction.
			if x[i][k]==1:
				for a,b in zip(R,L):
					nx=i+a
					ny=k+b
					#periodic boundary conditions
					if nx ==l_size:
						nx=0
					if nx==-1:
						nx=l_size-1
					if ny==l_size:
						ny=0
					if ny==-1:
						ny=l_size-1	
					if x[nx][ny]==0:
						n_1_0=n_1_0+1
					if x[nx][ny]==1:
						n_1_1=n_1_1+1
					if x[nx][ny]==2:
						n_1_2=n_1_2+1	
			if x[i][k]==2:
				for a,b in zip(R,L):
					nx=i+a
					ny=k+b
					if nx ==l_size:
						nx=0
					if nx==-1:
						nx=l_size-1
					if ny==l_size:
						ny=0
					if ny==-1:
						ny=l_size-1	
					if x[nx][ny]==0:
						n_2_0=n_2_0+1
					if x[nx][ny]==1:
						n_1_2=n_1_2+1
					if x[nx][ny]==2:
						n_2_2=n_2_2+1					
		
	#total the contribution from each type of interaction.	
	total_E_inter=n_1_2*eps_1_2+n_1_1*eps_1_1+n_2_2*eps_2_2
	
	return 	total_E_inter		

#function to find the total contribition due to vibrational energy by summing all terms in the vibrational energy matrix.
def calc_energy_vibr(E_vibr):
	tot_vibr_energy=sum(sum(E_vibr))
	
	return tot_vibr_energy
				

# function to find total energy of th wwhole system by summing up the contribution from each type of energy
def total_energy(x):
	tot_E=calc_enery_inter(x)+calc_energy_site(x)+calc_energy_vibr(E_vibr)
	return tot_E

#defining a function which moves the molecule on a chosen site in accordance with the rules of metropolis monte carlo, taking total energy into account	
def move_mol(cpx,cpy,x):
	if x[cpx][cpy]>0 and has_free(cpx,cpy)>0:
		orig_e=total_energy(x)
		n_list=nearest_neighbours(cpx,cpy)
		cs=randint(0,len(n_list[0])-1)
		x[n_list[0][cs]][n_list[1][cs]]=x[cpx][cpy]
		x[cpx][cpy]=0
		new_e=total_energy(x)
		if new_e-orig_e<0:
			x=x
		if 	new_e-orig_e>0:
			l=random.rand()
			acc_prob=min(1,np.exp(-beta*(new_e-orig_e)))
			if l>acc_prob:
				x=x
			else:
				x[cpx][cpy]=x[n_list[0][cs]][n_list[1][cs]]
				x[n_list[0][cs]][n_list[1][cs]]=0
					
	else:
		x=x		
		
	return x	
#define a function which changes state of a chosen site in accordance with metropolis Monte Carlo	
def change_state(cpx,cpy,x):
	if x[cpx][cpy]==1:
		orig_e=total_energy(x)
		x[cpx][cpy]=2
		new_e=total_energy(x)
		if new_e-orig_e<0:
			x=x
		if 	new_e-orig_e>0:
			l=random.rand()
			acc_prob=min(1,np.exp(-beta*(new_e-orig_e)))
			if l>acc_prob:
				x=x
			else:
				x[cpx][cpy]=1
	else:
		x=x
		
	return x				

#define a function which changes the vibrational energy of a chosen molecule if it is in state 2 in accordance with Metropolis Monte Carlo				
def change_vibr_energy(cpx,cpy,x,E_vibr):
	if x[cpx][cpy]==2:
		orig_vibr_e=E_vibr[cpx][cpy]
		orig_e=total_energy(x)
		y=random.rand()
		E_vibr[cpx][cpy]=E_vibr[cpx][cpy]+0.1*y*w
		new_e=total_energy(x)	
		if new_e-orig_e<0:
			x=x	
			E_vibr=E_vibr	
		if 	new_e-orig_e>0:
			l=random.rand()
			acc_prob=min(1,np.exp(-beta*(new_e-orig_e)))
			if l>acc_prob:
				x=x
				E_vibr=E_vibr
			else:
				E_vibr[cpx][cpy]=orig_vibr_e
	else:
		E_vibr=E_vibr
		x=x
	
	return E_vibr	
										



#setting the number of Monte Carlo Steps 		
Nm=10000

#choose a random point, and attempt one of the specified steps based on a probability.
for p in range(0,Nm):
	cpx=randint(0,l_size-1)
	cpy=randint(0,l_size-1)
	r=random.rand()
	if r<0.7:
		move_mol(cpx,cpy,x)
	if 0.7<r<0.85:
		change_state(cpx,cpy,x)
	if r>0.85: 
		change_vibr_energy(cpx,cpy,x,E_vibr)
		
plt.imshow(x)
plt.legend(x)
plt.show()

print(calc_energy_vibr(E_vibr))
print(calc_energy_site(x))
print(calc_enery_inter(x))
print(total_energy(x))
			
			

	


	
