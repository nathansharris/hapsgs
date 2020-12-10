#!/bin/python3
import numpy as np 
from itertools import combinations as comb
import vcf
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def getbreaks(frame,region):
    # Get the indices of the breaks. SGS breaks are places where opposite
    # heterozygotes occur. 
    allbreaks = np.array([i for i,j in enumerate(frame) if 0 in j and 2 in j])
    # Conver the indices to actual break positions. 
    locs = pos[allbreaks]
    # Get the breaks that occur inside of the given SGS region.
    inside  = locs[np.where( (locs >= region[0] ) & ( locs <= region[1]) )]
    # Get size of shared inside regions. 
    d = np.diff(inside)
    # Get the biggest one.
    largest = np.where(d == np.max(d))
    #specify a region that describes the largest shared region.
    out = np.array([inside[largest[0][0]],inside[largest[0][0]+1]])
    # We return the middle of the largest region, and the region.
    return np.array([np.mean(out),np.mean(out)]),out

def breakframe(frame,start,pos,direction):
    iter = np.arange(len(frame[0]))

    #ind = np.searchsorted(pos,start) #still having problems, and the mentality behind this line is 
    # the reason why. The Position of region are not always in pos because pos is only breaks.
    # We need to instead use the actual poition and not indices. this method assumes a match between
    # the index and the region position. This is not true. 
    # Compre all people pairwise. Look for places where genotypes sum to 2, but there are homozygotes. 
    #pairbreaks = np.array([ [np.where( (np.sum(frame[0:,[i,j]],axis=1) ==2) & ((frame[0:,i] == 0) | #pairbreaks is frame indices
    #    ( frame[0:,i] ==2 )) )[0] for i in iter ] for j in iter  ],dtype = object )
    pairbreaks = np.array([ [np.where( 
        (np.sum(frame[0:,[i,j]],axis=1) ==2) & (frame[0:,i] != 1 ) )[0] 
        for i in iter ] for j in iter  ],dtype = object )
    # dtype = object is new as far as I know. Numpy was upset without it. 

    # Now we are going to take the starting value and find its index in each set of pairwise breaks
    pairind = np.array( [[ np.searchsorted(pos[pairbreaks[i,j]], start) for i in iter ] # pairind is indices of a subset of pos. Using these
        for j in iter ] ) #indices to return to actual positions you need to use pos[pairbreaks[i,j]]
    
    nearest = np.array( [ [pos[pairbreaks[i,j]][pairind[i,j]] if i!=j else -1 for i in iter ] 
        for j in iter] )
    # nearest is the first break for each pair outside the SGS or interior region.
    if direction == "right":
        pass
    # or moving to the left
    else: # we need this block to adjust pairind when we want to move to the left, because searchsorted
        #will return the index to the right. 
        # Find breaks that do not land precisely on a region end.
        change = nearest != start
        # move their index one to the left. If this is being run on the left edge, this will
        # make nearest the first break to the left now. 
        pairind[change] -= 1
        nearest = np.array( [ [pos[pairbreaks[i,j]][pairind[i,j]-1] if i!=j else -1 for i in iter ] 
        for j in iter] )
    return nearest

# As far as I can tell this is called below once for an object that is not used. 
# I am leaving it in for now in case I remember why on earth it is here. 
def getorder(breaks): #breaks is distance to break from an edge.
    # Get all pairwise breaks locations
    vals = list(set(breaks[breaks>=0]))
    #we are going to iterate through them from smallest to largest
    vals.sort() 
    # we are going to use this and remove terms so we can maintain original order throughout.
    begin = np.arange(len(breaks))
    # Going to remove people one at a time.
    elimorder = []
    adj = []
    # iterate until two people left.
    while len(breaks) > 2:
        # Placeholder the length of the number of people. This will shrink each iteration.
        place = np.zeros(len(breaks))
        for val in vals: #walk through pairwise breaks starting with the nearest.
            # Find everyone who breaks with anyone at that particular break.
            place += np.sum(breaks == val, axis = 0)
            # If a column (person) has a break of this kind, add one to the same
            # index in the placeholder.
            # Once at least one person has four, stop. Why four...
            # Once at least one person has broken with everyone else, that person
            # and any other meeting the same condition. 
            if np.sum(place >= np.shape(breaks)[0]-1) > 0:
                break
                #The expression np.shape(breaks)[0]-1 used to be 4, but I think that is just when we
                # we using five people. Now it is number of remaining people minus one
        # get the indices of who is going to be removed based on their number of breaks.
        toremove = np.where(place == np.max(place))[0][0]
        elimorder.append(begin[toremove]) # add the index of the person(s) to the remove list.
        # Append the closest break to adj from the person(s) that were removed.
        adj.append(np.min(breaks[0:,toremove][breaks[0:,toremove] >0 ] ))
        # Take the removed person(s) out
        begin = np.delete(begin,toremove)
        # Removes the person from the breakframe, both as a row and a column.
        breaks = np.delete(np.delete(breaks,toremove,0),toremove,1)

    # Fill in the remaining two
    elimorder = elimorder+list(begin)
    adj = adj + [breaks[0,1],breaks[1,0]]

    # An array containing who broke first and where.
    s = np.array([elimorder,adj])
    # reorder to meet original ordering of people
    s = s[0:,np.argsort(s[0,0:])]
    # return its transpose
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
    
    # Print out the help message if no options are specified.
    if len(args) == 0:
        phelp()
    
    # This are the arguments that a user can use. 
    totalargs = ["--help","--interior","--skip","--pink","--regions","--file","--largest", "--solid"]
    
    # This will print out any arguments that are entered that are not in totalargs. 
    for i,arg in enumerate(args):
        if arg not in totalargs and args[i-1] not in totalargs:
            print("Argument " + str(arg) + " not understood.")
            print()
            print()
            phelp()
        else:
            continue

    # Allow the user to request the help message with a flag. 
    if "--help" in args:
        phelp()

    # Find the file flag and set the file variable as the following. 
    argin = np.where(args=="--file")[0]
    f = args[argin+1][0]

    # Do the same with the regions file. 
    argin = np.where(args=="--regions")[0]
    r = args[argin+1][0]

    # The option solid overrides the default changing color scheme of breaks. 
    try:
        argin = np.where(args=="--solid")[0]
        solcol = args[argin+1][0]
    except:
        solcol = "blue"

    # A default color scheme I am using. 
    colors = mcolors.TABLEAU_COLORS

    # Will use pay attention to breaks within the SGS region. 
    if "--interior" in args:
        interior = True
    else:
        interior = False

    # This is mostly for me. Setting the skip flag in conjunction with ipython's %run -i allows me to
    # skip the data import step o the run if the data is already imported. Makes for aster debugging. 
    if "--skip" in args:
        skip = True
    else:
        skip = False

    # The Julie option. Makes things pink. Very. 
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



    if skip == False: # Import data
        # This will generate a tabix file for the vcf if the user does not have one. 
        # Looking at you Julie.
        if not os.path.isfile(f+".tbi"):
            os.system("tabix "+f)


        # Read the regions file. 
        regions = np.genfromtxt(r, skip_header=1, dtype = int)
        # The regions file can contain one or multiple regions. Here we decide if
        # we need to iterate over multiple regions. 
        if len(regions.shape) == 1:
            regions = np.array([regions])[0:,1:]
        else:
            regions = regions[0:,1:]
        if len(np.unique(regions[0:,0])) > 1:
            print("Error. Multiple chrosomes specified in the regions file.") 
            print("Please be sure all regions are from the same chromosome and that these regions match your vcf.")
            sys.exit()
        else:
            chrom =  regions[0:,0][0]
        regions = regions[0:,1:]
            

        # Read in the vcf. One row per variant. 
        simple = np.array([record for record in vcf.Reader(filename=f).fetch(chrom,0)])
        
        # Retrieve snp positions.
        pos = np.array([rec.POS for rec in simple])

        # Convert vcf objects into 0,1,2 genotypes. One column per person, one
        # row per variant.
        simple = np.array([[ j.gt_type for j in i.samples ] for i in simple ])
        

        
  
    for ri,region in enumerate(regions):
        # Set a name for the ouput png.
        output=str(region[0])+"."+str(region[1])+"."+"png"

        # Plotting lines to show where the SGS region is. 
        plt.plot([regions[-1][0],regions[-1][0]], [-1,len(simple[0])], color = "gray",linestyle="--",label = "SGS")
        plt.plot([regions[-1][1],regions[-1][1]], [-1,len(simple[0])], color = "gray",linestyle="--")
        plt.ylim([-1,len(simple[0])])

        if interior == True:
            # Retrieve the inner region.
            region,subregion = getbreaks(simple,region)

            if largest == True:
                # Draws line around the largest interior region.
                plt.plot([subregion[0],subregion[0]], [-1,len(simple[0])], color = "black",linestyle="--",label = "Interior")
                plt.plot([subregion[1],subregion[1]], [-1,len(simple[0])], color = "black",linestyle="--")
                plt.ylim([-1,len(simple[0])])
        
        # Distance from first left break to region start
        leftbreaks = region[0] - breakframe(simple,region[0],pos,"left")
        np.fill_diagonal(leftbreaks,-1)
        # Distance from region end to first right break.
        rightbreaks = breakframe(simple,region[1],pos,"right") - region[1]
        np.fill_diagonal(rightbreaks,-1)

        # Get the entire chromosomal region for each person to last pairwise break on either side. 
        # Getorder returns distance to breaks. This is diabled for the time being. 
        #final = np.array ([region[0]-getorder(leftbreaks)[0:,1], getorder(rightbreaks)[0:,1]+region[1] ]).transpose()

            #print(np.sum(leftbreaks >= 0), np.sum(rightbreaks >= 0))

        # reorder by total. The idea here was to find the breaks on either size that was the furtherest away. We are looking
        # For the largest segment each person shares with at least one other person. 
        # new = np.mean(rightbreaks + leftbreaks, axis = 1).argsort()[::-1]

        # table of values we will use
        # right is now the location of breaks, not distance.
        right = (region[1] + rightbreaks)
        right = right.astype(float)
        np.fill_diagonal(right,np.nan)

        # Ditto for left. 
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

        np.savetxt("f"+output.strip("png")+"csv", ex, delimiter=",")
        np.savetxt("r"+output.strip("png")+"csv", right, delimiter=",")
        np.savetxt("l"+output.strip("png")+"csv", left, delimiter=",")


        
