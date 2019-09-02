import numpy as np
from matplotlib import pyplot as plt
import sys
import time 
from KeyPoller import KeyPoller


## in python 2, classes work differently (must inherite object in our case)
if sys.version_info[0] < 3:
    raise Exception("Python 3 is required, you are using ",
     sys.version_info[0])


class Gas:
    def __init__(self, amount_of_particles=50):
        self.__amount_of_particles = 0

        ## changeable parameters
        self._fs = 1                       # field size [0,_fs]x[0,_fs]
        self._vel = 0.015*self._fs         # update velocity 
        self._radius = 0.015*self._fs       # radius for seperation
        self._pause = 0.00001             # pause between steps must not be 0
        self._cos_angle = 0.0001             # cos(sight_angle)


        ## initialize array of particles
        self.__particles = np.empty(0, dtype="object")
        self.newParticles(amount_of_particles)

        plt.ion()
        self.__fig, self.__ax = plt.subplots()
        self.__sc = self.__ax.scatter(self.__position[0,:],self.__position[1,:], marker="o")
        plt.xlim(0,self._fs)
        plt.ylim(0,self._fs)

    def __updatePosDir(self):
        """
        update position and direction array for all particles
        self.__particles should be sorted before calling this function
        """
         ## array of x and y values are in the same order as self.__particles
        self.__position = np.empty([2,self.__amount_of_particles])
        self.__position[0,:] = np.array([ind.pos[0] for ind in self.__particles])
        self.__position[1,:] = np.array([ind.pos[1] for ind in self.__particles])

        ## direction array is in in the same order as self.__particles
        self.__direction = np.empty([2,self.__amount_of_particles])
        self.__direction[0,:] = np.array([ind.dir[0] for ind in self.__particles])
        self.__direction[1,:] = np.array([ind.dir[1] for ind in self.__particles])


    def update(self):
        """
        function: update position and direction of particles
        input:    none
        output:   none
        """
        ## ToDo: je nach iterationsrichtung aendert sich die Stromrichtung des Gases, 
        ## Stroeme nach oben/unten kommen nicht vor. Loesen, das hier ist nur provisorisch!!!
        marker = True
        if marker:
            for i in range(self.__amount_of_particles-1, -1, -1):
                self.__particles[i].calculateNewDirection(i)
            marker = False
        else:
            for i in range(0, self.__amount_of_particles):
                self.__particles[i].calculateNewDirection(i)
            marker = True
        for ind in self.__particles:
            ind.update()
        self.__updatePosDir()
        self.__particles = np.sort(self.__particles)

    
    def mapIndex(self, index):
        """
        input:    positive integer
        output:   index in range(0,self.__amount_of_particles) 
        """
        return int(index) % self.__amount_of_particles

    def distanceByIndex(self, ind1, ind2):
        """
        function: compute distance (manhattan metric) between two particles
        input:    index of particles in sorted list self.__particles
        output:   flaot, distance
        """
        ind1 = self.mapIndex(ind1)
        ind2 = self.mapIndex(ind2)

        dist_x = abs(self.__particles[ind1].pos[0]-self.__particles[ind2].pos[0])
        if (dist_x>self._fs/2):
            dist_x = self._fs-dist_x
        dist_y = abs(self.__particles[ind1].pos[1]-self.__particles[ind2].pos[1])
        if (dist_y>self._fs/2):
            dist_y = self._fs-dist_y
        return dist_x+dist_y
                
    
    def xDistByIndex(self, ind1, ind2):
        """
        function: test if index is in range of (0, self.__amount_of_particles)
        input:    list or set of indices of particles in sorted array self.__particles
        output:   bool, True if in range, False if not
        """
        ind1 = self.mapIndex(ind1)
        ind2 = self.mapIndex(ind2)
        dist_x = abs(self.__particles[ind1].pos[0]-self.__particles[ind2].pos[0])
        if (dist_x>self._fs/2):
            dist_x = self._fs-dist_x
        return dist_x

    def inRange(self, index):
        """
        function: test if index is in range of (0, self.__amount_of_particles)
        input:    list or set of indices of particles in sorted array self.__particles
        output:   bool, True if in range, False if not
        """
        try: 
            i = int(index)
            return (i < self.__amount_of_particles and i >= 0)
        except:
            return all(x in range(0,self.__amount_of_particles) for x in index)


    def getParticleByIndex(self, index):
        """
        function: get func
        input:    index of particle in sorted array self.__particles
        output:   none
        """
        assert self.inRange(index), "index not in range"
        return self.__particles[int(index)]


 
    def angleByIndex(self, ind_self, ind_other):
        """
        function: test if bird "other" is in sight range of bird "self"
        input:    indices of individuals in sorted array self.__individuals
        output:   bool, True if in range, False if not
        """
        ret_val = False
        ind_self = self.mapIndex(ind_self)
        ind_other = self.mapIndex(ind_other)

        a = self.getParticleByIndex(ind_self).dir
        b = self.getParticleByIndex(ind_other).pos-self.getParticleByIndex(ind_self).pos
        try: 
            ret_val = np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)) > self._cos_angle
        except:
            ## case if one of the vecor norms is zero
            ret_val = True
        return ret_val



    class Particle:
        """ 
        nested class containing particles (agents) and their properties
        """
        def __init__(self, gas):
            self.g = gas
            ## define position pos[0]..x, pos[1]..y
            self.pos = np.array([float(np.random.random()*self.g._fs), float(np.random.random()*self.g._fs)])
            ## define current direction
            self.dir = np.array([float(np.random.uniform(-1,1,1)), float(np.random.uniform(-1,1,1))])
            self.partner = None
            self.new_dir_calculated = False

        ## functions to compare objects
        def __eq__(self, other):
            return eq(self.pos[0], other.pos[0])
        def __ne__(self, other):
            return ne(self.pos[0], other.pos[0])
        def __lt__(self, other):
            return self.pos[0] < other.pos[0]
        def __le__(self, other):
            return self.pos[0] <= other.pos[0]
        def __gt__(self, other):
            return self.pos[0] > other.pos[0]
        def __ge__(self, other):
            return self.pos[0] >= other.pos[0]
        
        def update(self):
            """
            function: update direction and position of particle
            input:    none
            output:   none
            """
            self.dir = self.new_dir
            self.pos[0] = (self.pos+self.new_dir*self.g._vel)[0] % self.g._fs
            self.pos[1] = (self.pos+self.new_dir*self.g._vel)[1] % self.g._fs
            self.partner = None
            self.new_dir_calculated = False

        def calculateNewDirection(self, index):
            """
            function: calculate new direction for this particle and store it in self.new_dir
            input:    index of this particle in sorted array gas.__particle
            return:   none
            """

            if self.partner == None:
                i = index

                ## find particles in environment
                ## inc = -1: go to the left in s.particles list, inc = 1: go to the right
                increment = np.random.choice([-1,1])
                for inc in [increment, -increment]:
                    i = index+inc
                    while self.g.mapIndex(i) != index and self.g.xDistByIndex(index, i) <= self.g._radius:
                        if self.g.distanceByIndex(index, i) <= self.g._radius and self.g.getParticleByIndex(self.g.mapIndex(i)).partner == None and self.g.angleByIndex(index, i):
                            self.g.getParticleByIndex(self.g.mapIndex(i)).partner = self.g.mapIndex(index)
                            self.partner = self.g.mapIndex(i)
                            break
                        i+=inc
    
            if (not self.new_dir_calculated) and self.partner != None:
                if self > self.g.getParticleByIndex(self.partner):
                    greater = self
                    smaller = self.g.getParticleByIndex(self.partner)          
                else:
                    smaller = self
                    greater = self.g.getParticleByIndex(self.partner)                      

                n = greater.pos - smaller.pos  ## stossgerade
                #input("""stossgerade: {}, greater.pos: {}, grater.dir: {} 
                #smaller.pos: {}, smaller.dir: {}""".format(n, greater.pos, greater.dir, smaller.pos, smaller.dir))
                cos_phi = np.dot([1,0],n)/np.linalg.norm(n)
                sin_phi = np.sqrt(1-cos_phi**2)
                if greater.pos[1] >= smaller.pos[1]:
                    dir_of_rotation = -1 ## counter clockwise
                else:
                    dir_of_rotation = 1 ## clockwise
                
                rot_mat = np.array([[cos_phi, -dir_of_rotation*sin_phi], [dir_of_rotation*sin_phi, cos_phi]]) ## rotation matrix
                
                ## simulate collision with same weights 
                rot_gr = np.matmul(rot_mat, greater.dir)
                rot_sm = np.matmul(rot_mat, smaller.dir)
                #input("""NEU stossgerade: {}, greater.pos: {}, grater.dir: {} 
                #smaller.pos: {}, smaller.dir: {}""".format(np.matmul(rot_mat, n), np.matmul(rot_mat, greater.pos), rot_gr, 
                #np.matmul(rot_mat, smaller.pos), rot_sm))

                tmp = rot_gr[0]
                rot_gr[0] = rot_sm[0]
                rot_sm[0] = tmp

                rot_mat_inv = [[cos_phi, dir_of_rotation*sin_phi], [-dir_of_rotation*sin_phi, cos_phi]] ## inverse rotation matrix
                greater.new_dir = np.matmul(rot_mat_inv, rot_gr)
                smaller.new_dir = np.matmul(rot_mat_inv, rot_sm)

                greater.new_dir_calculated = True
                smaller.new_dir_calculated = True
            else:
                self.new_dir = self.dir
                self.new_dir_calculated = True


    def newParticles(self, amount_of_particles=1):
        """
        function: create amount_of_particles new particles and add them to self.__particles list
        input:    amount_of_particles, integer, default 1
        output:   none
        """
        assert int(amount_of_particles)>0, "must be positive integer"
        x_lex = [self.Particle(self) for i in range(int(amount_of_particles))]
        self.__amount_of_particles += int(amount_of_particles)
        self.__particles = np.sort(np.append(self.__particles, np.array(x_lex, dtype="object")))
        self.__updatePosDir()
        print(amount_of_particles, " new particles created")

    
    def delParticles(self, amount_of_particles=1):
        """
        function: delete amount_of_particles particles from self.__particles list
        input:    amount_of_particles, integer, default 1
        output:   none
        """
        assert int(amount_of_particles)>0, "must be positive integer"
        amount_of_particles = int(amount_of_particles)
        if amount_of_particles >= self.__amount_of_particles:
            print("delete all birds except one instead of {}".format(amount_of_particles))
            amount_of_particles = self.__amount_of_particles-1
        to_del = np.random.choice(self.__amount_of_particles, size=amount_of_particles, replace=False)
        # ToDo: make sure destructor is called (?)
        rem = np.delete(self.__particles, to_del)
        self.__particles = np.sort(rem)
        self.__amount_of_particles -= amount_of_particles
        self.__updatePosDir()  
        print(amount_of_particles, " particles deleted")


    def plot(self):
        """
        function: plot gas in scatter plot (each particle is a dot)
        input:    none
        output:   none
        """
        self.__sc.set_offsets(np.c_[self.__position[0,:],self.__position[1,:]])
        self.__fig.canvas.draw_idle()
        if self._pause >0: plt.pause(self._pause)
    
    def changeParameters(self):
        """
        function: change changeable parameters via comand line inputs
        input:    none
        output:   none
        """
        print("""If you have finished changing parameters press [q] + [ENTER]
        Current Parameters:
        [1] field size: {}
        [2] velocity: {}
        [3] seperation radius (as factor of field size): {}*field_size
        [4] alignment radius (as factor of field size): {}*field_size
        [5] cohesion radius (as factor of field size): {}*field_size
        [6] seperation weight: {}
        [7] alignment weight: {}
        [8] cohesion weight: {}
        [9] pause between steps: {}
        [10] cos of range of vision: {}"""
        .format(self._fs, self._vel, self._radius/self._fs, self._radius_ali/self._fs, 
        self._radius_co/self._fs, self._weight_sep, self._weight_ali, self._weight_co, self._pause, self._cos_angle))

        while True:
            no = input("type number of parameter you want to change + [ENTER] or [q] + [ENTER] to escape \n")
            if no == "q": break
            if no == "1":
                value = input("type value + [ENTER] \n")
                try: 
                    if int(value)>0:
                        self._fs = int(value)
                        plt.xlim(0,self._fs)
                        plt.ylim(0,self._fs)
                    else:
                        print("value should be positive integer")
                except:
                    print("value should be positive integer")
            elif no == "2":
                value = input("type value + [ENTER] \n")
                try: 
                    if float(value)>0:
                        self._vel = float(value)
                    else:
                        print("value should be positive float")
                except:
                    print("value should be positive float")
            elif no == "3":
                value = input("type value + [ENTER] \n")
                try: 
                    if float(value)>0 and float(value)<=1 and float(value)<=self._radius_ali:
                        self._radius = float(value)*self._fs
                    else:
                        print("value should be positive float that is smaller than 1 and smaller than alignment radius")
                except:
                    print("value should be positive float that is smaller than 1 and smaller than alignment radius")
            elif no == "4":
                value = input("type value + [ENTER] \n")
                try: 
                    if float(value)>0 and float(value)<=1 and float(value)<=self._radius_co and float(value)>=self._radius:
                        self._radius_ali = float(value)*self._fs
                    else:
                        print("value should be positive float that is smaller than 1 and between seperation and cohesion radius")
                except:
                    print("value should be positive float that is smaller than 1 and between seperation and cohesion radius")
            elif no == "5":
                value = input("type value + [ENTER] \n")
                try: 
                    if float(value)>0 and float(value)<=1 and float(value)>=self._radius_ali:
                        self._radius_co = float(value)*self._fs
                    else:
                        print("value should be positive float that is smaller than 1 and geater than seperation radius")
                except:
                    print("value should be positive float that is smaller than 1 and greater than seperation radius")
            elif no == "6":
                value = input("type value + [ENTER] \n")
                try:
                    self._weight_sep = float(value)
                except:
                    print("value should be float")
            elif no == "7":
                value = input("type value + [ENTER] \n")
                try:
                    self._weight_ali = float(value)
                except:
                    print("value should be float")
            elif no == "8":
                value = input("type value + [ENTER] \n")
                try:
                    self._weight_co = float(value)
                except:
                    print("value should be float")
            elif no == "9":
                value = input("type value + [ENTER] \n")
                try:
                    if float(value) > 0: 
                        self._pause = float(value)
                    else:
                        print("value should be positive float")
                except:
                    print("value should be positive float")
            elif no == "10":
                value = input("type value + [ENTER] \n")
                try:
                    if float(value) > -1 and float(value) < 1: 
                        self._pause = float(value)
                    else:
                        print("value should be float between -1 and 1")
                except:
                    print("value should be float between -1 and 1")
            else:
                print("invalid input")

amount = input("How many particles should the gas contain? input [amount of particles] + [ENTER] \n")

## create gas
gas = Gas(amount_of_particles=int(amount))
with KeyPoller() as keyPoller:
    print("""if you want to end the program press [q]
    if you want to change the gas parameters press [p] 
    if you want to add particles press [a] 
    if you want to delete particles press [d] 
    (there are still some bugs, I'll try to fix, sorry for that)""")
    
    while True:
        tmp = keyPoller.poll()
        if not tmp is None:
            if tmp == "q":
                break
            elif tmp== "p":
                gas.changeParameters()
            elif tmp == "a":
                value = input("input [amount to add] + [ENTER]\n")
                try: 
                    if int(value)>0: 
                        gas.newParticles(int(value))
                    else: 
                        print("input should be positive integer")
                except:
                    print("input should be positive integer")
            elif tmp=="d":
                value = input("input [amount to delete] + [ENTER]\n")
                try: 
                    if int(value)>0: 
                        gas.delparticles(int(value))
                    else: 
                        print("input should be positive integer")
                except:
                    print("input should be positive integer")
            else:
                print("invalid input")
                
        gas.update()
        gas.plot()







