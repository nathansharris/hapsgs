#!/bin/python3
import numpy as np 
from itertools import combinations as comb
import vcf
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def getbreaks(frame,region):
    allbreaks = np.array([i for i,j in enumerate(frame) if 0 in j and 2 in j])
    locs = pos[allbreaks]
    inside  = locs[np.where( (locs >= region[0] ) & ( locs <= region[1]) )]
    d = np.diff(inside)
    largest = np.where(d == np.max(d))
    out = np.array([inside[largest[0][0]],inside[largest[0][0]+1]])
    return np.array([np.mean(out),np.mean(out)]),out

def breakframe(frame,start,pos,dir):
    iter = np.arange(len(frame[0]))

    #ind = np.searchsorted(pos,start) #still having problems, and the mentality behind this line is 
    # the reason why. The Position of region are not always in pos because pos is only breaks.
    # We need to instead use the actual poition and not indices. this method assumes a match between
    # the index and the region position. This is not true. 
    pairbreaks = np.array([ [np.where( (np.sum(frame[0:,[i,j]],axis=1) ==2) & ((frame[0:,i] == 0) | #pairbreaks is frame indices
        ( frame[0:,i] ==2 )) )[0] for i in iter ] for j in iter  ] )

    pairind = np.array( [[ np.searchsorted(pos[pairbreaks[i,j]], start) for i in iter ] # pairind is indices of pairind
        for j in iter ] )
        
    nearest = np.array( [ [pos[pairbreaks[i,j]][pairind[i,j]] if i!=j else -1 for i in iter ] 
        for j in iter] )
    
    if dir == "right":
        pass
    if dir == "left": # we need this block to adjust pairind when we want to move to the left, because searchsorted
        change = nearest != start #will return the index to the right. 
        pairind[change] -= 1
        nearest = np.array( [ [pos[pairbreaks[i,j]][pairind[i,j]-1] if i!=j else -1 for i in iter ] 
        for j in iter] )
    return nearest

def getorder(breaks):
    vals = list(set(breaks[breaks>=0]))
    vals.sort()
    begin = np.arange(len(breaks))
    elimorder = []
    adj = []
    while len(breaks) > 2:
        place = np.zeros(len(breaks))
        for val in vals: 
            for throw in np.concatenate(np.where(breaks == val)):
                place[throw]+=1
            if np.sum(place >= 4) > 0:
                break
        toremove = np.where(place == np.max(place))[0][0]
        elimorder.append(begin[toremove])
        adj.append(np.min(breaks[0:,toremove][breaks[0:,toremove] >0 ] ))
        begin = np.delete(begin,toremove)
        breaks = np.delete(np.delete(breaks,toremove,0),toremove,1)

    elimorder = elimorder+list(begin)
    adj = adj + [breaks[0,1],breaks[0,1]]

    s = np.array([elimorder,adj])
    s = s[0:,np.argsort(s[0,0:])]

    return s.transpose()

def phelp():
    print("This script outputs a visual representation of haplotype sharing around a given SGS region.")
    print()
    print("Required arguments:")
    print("    --file <vcf> \tVCF file containing the individuals you wish to display.")
    print("    --regions <txt> \tA plain text file with the regions you wish to focus on")
    print()
    print("Optional arguments:")
    print("    --interior \t\tSetting this flag will start the haplotype sharing at the" )
    print("    \t\t\tlargest shared region inside the SGS region.")
    print("    --pink\t\tUses a single color, right now pink, with a decreasing")
    print("\t\t\talpha instead of many")
    print("    --largest\t\tShow the largest interior shared region. --interior")
    print("\t\t\tflag must be set.")

    sys.exit()

if __name__ == "__main__":
    args = np.array(sys.argv[1:])
    if len(args) == 0:
        phelp()
    
    totalargs = ["--help","--interior","--skip","--pink","--regions","--file","--largest", "--solid"]
    
    for i,arg in enumerate(args):
        if arg not in totalargs and args[i-1] not in totalargs:
            print("Argument " + str(arg) + " not understood.")
            print()
            print()
            phelp()
        else:
            continue


    if "--help" in args:
        phelp()


    argin = np.where(args=="--file")[0]
    f = args[argin+1][0]
    argin = np.where(args=="--regions")[0]
    r = args[argin+1][0]
    try:
        argin = np.where(args=="--solid")[0]
        solcol = args[argin+1][0]
    except:
        solcol = "blue"

    colors = mcolors.TABLEAU_COLORS

    if "--interior" in args:
        interior = True
    else:
        interior = False
    if "--skip" in args:
        skip = True
    else:
        skip = False
    if "--pink" in args:
        coloption = "alpha"
    elif "--solid" in args:
        coloption = "solid"
    else:
        coloption = "no"
    if "--largest" in args:
        largest = True
    else:
        largest = False



    if skip == False: #optionally skip data import
        if not os.path.isfile(f+".tbi"):
            os.system("tabix "+f)
        simple = np.array([record for record in vcf.Reader(filename=f).fetch(2,0)])
        pos = np.array([rec.POS for rec in simple])
        simple = np.array([[ j.gt_type for j in i.samples ] for i in simple ])
        
        regions = np.genfromtxt(r, skip_header=1, dtype = int)
        
        if len(regions.shape) ==1:
            regions = np.array([regions])[0:,1:]
        else:
            regions = regions[0:,1:]
        
  
    for ri,region in enumerate(regions):
        output=str(region[0])+"."+str(region[1])+"."+"png"

        plt.plot([regions[-1][0],regions[-1][0]], [-1,len(simple[0])], color = "gray",linestyle="--",label = "SGS")
        plt.plot([regions[-1][1],regions[-1][1]], [-1,len(simple[0])], color = "gray",linestyle="--")
        plt.ylim([-1,len(simple[0])])

        if interior == True:
            region,subregion = getbreaks(simple,region)

            if largest == True:
                plt.plot([subregion[0],subregion[0]], [-1,len(simple[0])], color = "black",linestyle="--",label = "Interior")
                plt.plot([subregion[1],subregion[1]], [-1,len(simple[0])], color = "black",linestyle="--")
                plt.ylim([-1,len(simple[0])])
        
        leftbreaks = region[0] - breakframe(simple,region[0],pos,"left")
        np.fill_diagonal(leftbreaks,-1)
        rightbreaks = breakframe(simple,region[1],pos,"right") - region[1]
        np.fill_diagonal(rightbreaks,-1)

        final = np.array ([region[0]-getorder(leftbreaks)[0:,1], getorder(rightbreaks)[0:,1]+region[1] ]).transpose()

            #print(np.sum(leftbreaks >= 0), np.sum(rightbreaks >= 0))

        #reorder by total 
        new = np.mean(rightbreaks + leftbreaks,axis = 1).argsort()[::-1]
        rightbreaks = rightbreaks[new]
        leftbreaks = leftbreaks[new]

        # table of values we will use
        right = (rightbreaks + region[1])
        right = right.astype(float)
        np.fill_diagonal(right,np.nan)

        left = (region[0] - leftbreaks)
        left = left.astype(float)
        np.fill_diagonal(left,np.nan)

        farright = np.nanmax(right, axis = 1)
        farleft = np.nanmax(left, axis = 1)
        ex = np.array([farleft, farright]).T


        #graphing with solid color
        if "--solid" in args:
            for i,row in enumerate(ex):
                plt.plot([row[0],row[1]], [i,i], color = solcol,linewidth=9,solid_capstyle="butt")

        else:
            #Graph with changing colors
            for i,person in enumerate(rightbreaks):
                sharing = len(simple[0])-1
                a = 1
                if coloption == "alpha":
                    ourcolor = "magenta"
                elif coloption == "solid":
                    ourcolor = "blue"
                else:
                    ourcolor = list(colors)[sharing]
                plt.plot([region[0],region[1]],[i,i],color = ourcolor,linewidth = 9,solid_capstyle="butt")
                start = region[1]
                for j,pair in enumerate(sorted(list(set(person[person>= 0])))):
                    numbreaking = np.sum(person == pair)
                    if coloption == "alpha":
                        ourcolor = "magenta"
                        alpha = 1 - a*(1/(len(simple[0])))
                    elif coloption == "solid":
                        ourcolor = "blue"
                        alpha = 1
                    else:
                        ourcolor = list(colors)[sharing]
                        alpha = 1
                    plt.plot([start,region[1]+pair], [i,i], color = ourcolor,linewidth=9,label = str(sharing),solid_capstyle="butt",alpha =alpha)
                    sharing -= numbreaking
                    a+= numbreaking
                    start = region[1]+pair

            for i,person in enumerate(leftbreaks):
                sharing = len(simple[0])-1
                a=1
                start = region[0]
                for j,pair in enumerate(sorted(list(set(person[person>= 0])))):
                    numbreaking = np.sum(person == pair)
                    if coloption == "alpha":
                        ourcolor = "magenta"
                        alpha = 1 - a*(1/(len(simple[0])))
                        print(alpha)
                    else:
                        ourcolor = list(colors)[sharing]
                        alpha = 1
                    plt.plot([region[0]-pair,start], [i,i], color = ourcolor,linewidth=9,label = str(sharing),solid_capstyle="butt",alpha = alpha  )
                    sharing -= numbreaking
                    a += numbreaking 
                    start = region[0]-pair

        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        
        order = np.argsort(list(by_label.keys()))
        keys = np.array(list(by_label.keys()))[order]
        vals = np.array(list(by_label.values()))[order] 
        plt.legend(vals, keys,title="Sharing")
        plt.yticks([])
        #plt.xlabel("Position (bp)")
        plt.savefig(output,dpi = 400)
        plt.show()


        
