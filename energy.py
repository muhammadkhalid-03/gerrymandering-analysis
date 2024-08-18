# energy.py
# Eric A. Autry
# CSC 395 Spring 2024
# 05/09/24

import math
import numpy as np

"""
weightedMu: calculates a list of energies and takes a weighted sum.

energyLst: the energy info as a list of (weight, energyFunc) pairs
plan: the initial districting plan matrix
border: the initial list of current border nodes (row,col) pairs

returns: the weighted sum of each energy value
         sum( weight*energyFunc(plan,border) for each energy )
"""
def weightedMu(energyLst, plan, border):
  energy = 0
  for weight, energyFunc in energyLst:
    energy += weight*energyFunc(plan, border)
  return math.exp(-energy)

"""
getHVmatrices: get the H and V matrices for the (a)symmetry energies.
               Note: better to pre-compute and reuse.

N: size of the districting plan (NxN)

returns: H, V matrices of the form

H = [ 1  1  1  1  1 
      1  1  1  1  1
      0  0  0  0  0 
     -1 -1 -1 -1 -1
     -1 -1 -1 -1 -1 ]

V = H.T
"""
def getHVmatrices(N):
  halfN = math.floor(N/2)
  V = np.ones((N,N))
  V *= -1
  V[:, :halfN] = 1
  if (N%2 == 1):
    V[:, halfN] = 0
  H = V.copy()
  H = H.T
  return H, V

# Below are the four (a)symmetry energy functions.
# horz_sym: measure for how symmetric the plan is about the horizontal axis
# vert_sym: measure for how symmetric the plan is about the vertical axis
# horz_asym: measure for how asymmetric the plan is about the horizontal axis
# vert_asym: measure for how asymmetric the plan is about the vertical axis
#
# Inputs: current plan X, horizontal matrix H, vertical matrix V.
# 
# H = [ 1  1  1  1  1 
#       1  1  1  1  1
#       0  0  0  0  0 
#      -1 -1 -1 -1 -1
#      -1 -1 -1 -1 -1 ]
# 
# V = H.T
#
# A note about these energies: recall that we will use a mu of the form 
# e^-energy, so we will tend to see smaller energy plans. These energies are
# all practically calculated with something like sum(X*V). Since X is a metrix
# of 1s and 0s, and V (or H) is a matrix of -1s, 0s, 1s, then their product
# will be minimized when the 1s of X align with the -1s of H:
# horizontal asymmetry.
#
# However, when we take the absolute value, we will want the 1s and -1s to 
# cancel. This will mean that the H matrix will actually prioritize balancing
# the 1s of X between the top and bottom of the domain: horizontal symmetry.
#
# tl, dr
# The asym energies will tend to align the 1s in the plan X to the -1s in
# whichever matrix V or H is involved.
# The sym energies will try to balance 1s and -1s in H or V.
def horz_sym(X, H, V):
  return abs(np.sum(X*H))

def vert_sym(X, H, V):
  return abs(np.sum(X*V))

def horz_asym(X, H, V):
  return np.sum(X*H)

def vert_asym(X, H, V):
  return np.sum(X*V)

"""
isoperimetric: compute the isoperimetric score used as a measure of 
               compactness. This score is calculated as:
               perimeter^2 / area

               When used as an energy, we will tend to minimize the perimeter 
               with respect to the area. The actual minimizing shape is a 
               circle, so this measure will try to make shapes that look...
               somewhat circular ('blobs'). It prefers smoother boundaries.

               To calculate the these geometric quantities, we will set our 
               square precincts to be 1x1 unit squares with area 1 and perim 4.

               The area is then the number of 0s (or 1s) in the plan.

               The perimeter is the number of 0s (or 1s) on the border plus
               the number of 0s (or 1s) along the edges of the overall domain
               (i.e., in the top/bottom row or left/right column). These count
               as boundary nodes as well, so contribute to overall perimeter 
               of the districts.

plan: the initial districting plan matrix
N: the size of the plan (NxN)
border: the initial list of current border nodes (row,col) pairs

returns: the average isoperimetric ratio of the two districts in the plan
"""
def isoperimetric(plan, N, border):
  # List of vals (district labels) to avg over.
  val_lst = [0,1]

  # List to store the value per district.
  isoperimetric_lst = [0,0]

  # For each district (val).
  for ii in range(len(val_lst)):
    val = val_lst[ii]

    # The area of the dist is the count of the labelled nodes.
    area = sum(sum(plan==val))

    # The perimeter!
    perim = 0.0

    # Non-corners (have 1 perim)
    perim += sum(plan[0,1:N-2]==val)
    perim += sum(plan[N-1,1:N-2]==val)
    perim += sum(plan[1:N-2,0]==val)
    perim += sum(plan[1:N-2,N-1]==val)

    # Corners (have 2 perim)
    perim += 2*(plan[0,0]==val)
    perim += 2*(plan[N-1,0]==val)
    perim += 2*(plan[0,N-1]==val)
    perim += 2*(plan[N-1,N-1]==val)

    # Count the boundary node's neighbors (+1 each).
    for (br, bc) in border:
      if plan[br,bc] == val:
        if br > 0:
          if plan[br-1,bc] != val:
            perim += 1
        if br < N-1:
          if plan[br+1,bc] != val:
            perim += 1
        if bc > 0:
          if plan[br,bc-1] != val:
            perim += 1
        if bc < N-1:
          if plan[br,bc+1] != val:
            perim += 1

    # The score.
    isoperimetric_lst[ii] = (perim**2) / area

  # Avg over the districts (vals).
  return sum(isoperimetric_lst)/len(val_lst)
