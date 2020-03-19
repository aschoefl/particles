#!python3
import sys, random
import pymunk
import numpy as np 
from scipy import constants

import pygame
from pygame.locals import *
import pymunk.pygame_util

## graphic parameters
GRAPHICS = True # True - graphic output (slow); False - no graphic output (faster)
DT = 0.01        # timestap parameter for graphics

## physical parameters
T = 410             # [K]
k = 1.38064852    # actually 1.38064852*10^-23 [J/K]
AMOUNT_OF_NO = 500  
AMOUNT_OF_O3 = 500
RADIUS_OF_PARTICLES = 1   # actually *10^-10 m
MASS_OF_NO = (2.3+2.67)/10**3  # actually *10^-26 kg
MASS_OF_O3 = 2.67*3/10**3     # actually *10^-26 kg
VELOCITY_MEAN_NO = (k*T*constants.pi/MASS_OF_NO/2)**(1/2)  
VELOCITY_STD_NO = (k*T/MASS_OF_NO*(1-constants.pi/4))**(1/2)  
VELOCITY_MEAN_O3 = (k*T*constants.pi/MASS_OF_O3/2)**(1/2)
VELOCITY_STD_O3 = (k*T/MASS_OF_O3*(1-constants.pi/4))**(1/2)  
THRESHOLD_ENERGY = 1.93*10**3  # actually 1,93*10^-20 J
REACTION_ENERGY = 3.093*10**2  # actually 3,093*10^-19 J

## model parameters
MAX_TIME = 500      # time per round
MAX_ROUNDS = 300    # amount of runs
SPACE_SIZE = 1000  
ELASTICITY = 0.999  # should not be >= 1 according to docu of pymunk
LINES_RADIUS = 4 

## set random seed
random.seed = 42

## initiate global counter
amount_of_reactions = 0

##############################
## 
##############################
def addParticles(space, amount_list, molecule_type_list, existing_particles = []):
    """ 
    function that adds particles for an existing space

    Parameters:
    -----------
    space: space object as returned by pymunk.Space()
    amount_list: list of int
        containing the amount of particles to be produced for each
        type of molecule in molecule_type_list
    molecule_type_list: list of str
        containig types of particles to be produced (allowed: "NO" and "O3")
        in the same order as in amount_list
    existing_particles: list of objects
        list of all particles already existing in the space

    Returns:
    --------
    list of all particles in the space. contains objects as
    returned by pymunk.Body(.)
    """
    try:
        amount_list = list(amount_list)
    except:
        tmp = amount_list
        amount_list = []
        amount_list.append[tmp]
    try:
        molecule_type_list = list(molecule_type_list)
    except:
        tmp = molecule_type_list
        molecule_type_list = []
        molecule_type_list.append[molecule_type_list]

    assert len(molecule_type_list)==len(amount_list)
    particles = existing_particles
    
    for i in range(0, len(amount_list)):
        molecule_type = molecule_type_list[i]
        amount = amount_list[i]
        """Add a ball to the given space at a random position"""
        if molecule_type == "NO":
            mass = MASS_OF_NO
            vel_mean = VELOCITY_MEAN_NO
            vel_std = VELOCITY_STD_NO
        elif molecule_type == "O3":
            mass = MASS_OF_O3
            vel_mean = VELOCITY_MEAN_O3
            vel_std = VELOCITY_STD_O3
        else:
            assert("wrong input for molecule_type")
        radius = RADIUS_OF_PARTICLES

        for j in range(0, amount):
            inertia = pymunk.moment_for_circle(mass, 0, radius, (0,0))
            body = pymunk.Body(mass, inertia)
            ## set initial impuls and coordinates of particles
            ## gaussian distribution of velocity
            p_abs = random.gauss(vel_mean, vel_std)*mass
            ## uniform distribution of coordinates
            x_dir = random.random()*random.choice([-1,1])
            y_dir = (1-x_dir**2)**(1/2)*random.choice([-1,1])
            body.apply_impulse_at_local_point((x_dir*p_abs, y_dir*p_abs))
            x = random.uniform(radius/2+LINES_RADIUS/2+2,SPACE_SIZE-radius/2-LINES_RADIUS/2-2)
            y = random.uniform(radius/2+LINES_RADIUS/2+2,SPACE_SIZE-radius/2-LINES_RADIUS/2-2)
            body.position = x, y
            shape = pymunk.Circle(body, radius, (0,0))
            shape.elasticity = ELASTICITY
            shape.friction = 0
            ## definition of collision_type: 0-lines, 1-particle
            shape.collision_type = 1
            shape.molecule_type = molecule_type
            space.add(body, shape)
            try:
                particles.append(body)
            except:
                "Wrong type for existing_particles in addParticles"
    return(particles)


def addLines(space):
    """ 
    function that draws a quadratic outline of size SPACE_SIZExSPACE_SIZE
    
    Parameters:
    -----------
    space: space object as returned by pymunk.Space()

    """

    body = pymunk.Body(body_type = pymunk.Body.STATIC)
    body.position = (0,0)
    body.elasticity = ELASTICITY

    lines = []

    lines.append(pymunk.Segment(body, (0, 0), (SPACE_SIZE, 0.0), LINES_RADIUS))
    lines.append(pymunk.Segment(body, (0, 0), (0.0, SPACE_SIZE), LINES_RADIUS))
    lines.append(pymunk.Segment(body, (SPACE_SIZE, 0), (SPACE_SIZE, SPACE_SIZE), LINES_RADIUS))
    lines.append(pymunk.Segment(body, (0, SPACE_SIZE), (SPACE_SIZE, SPACE_SIZE), LINES_RADIUS))


    for l in lines:
        l.elasticity = ELASTICITY
        l.friction = 0
        ## definition of collision_type: 0-lines, 1-particle
        l.collision_type = 0
    space.add(lines[0], lines[1], lines[2], lines[3], body)

def reaction(v1, v2, m1, m2):
    """ 
    determines if relative kinetic energy is high enough
    for a reaction to happen
    
    Parameters:
    -----------
    v1, v1: int
        velocity of colliding objects
    m1, m2: int
        mass of colliding objects

    Returns:
    --------
    True if collision happens, False else

    """
    if m1*m2*np.dot(np.array(v1-v2),np.array(v1-v2))/(2*(m1+m2)) >= THRESHOLD_ENERGY:
        return True
    else:
        return False

def set_new_velocites(p1, p2):
    """ 
    simulates a exothermic reaction by altering the particle's velocity

    Parameters:
    -----------
    p1, p2: objects as returned by pymunk.Body(.)
        colliding objects
    """
    # divide ekin by 2 beacause apparently:
    # s[j].body.kinetic_energy === np.linalg.norm(s[j].body.velocity)**2*s[j].body.mass
    alpha = 1 + REACTION_ENERGY/(p1.body.kinetic_energy/2+p2.body.kinetic_energy/2)
    p1.body.velocity = p1.body.velocity*alpha
    p2.body.velocity = p2.body.velocity*alpha
    
def collision_handler(arbiter, space, data):
    """ 
    pre solve collision handler (see pymunk documentation for details)
    determines if a reaction between the colliding objects is happening

    Parameters:
    -----------
    arbiter: pymunk arbiter object
    space: space object as returned by pymunk.Space()
    """
    global amount_of_reactions
    s = arbiter.shapes
    #n = np.array(arbiter.contact_point_set.normal)
    for j in range(0, len(s)):
        if s[j].molecule_type in ["NO2, O2"]:
            pass
        elif s[j].molecule_type == "NO":
            for i in range(j+1, len(s)):
                if s[i].molecule_type == "O3" and reaction(s[j].body.velocity, s[i].body.velocity,
                            s[j].body.mass, s[i].body.mass):
                    ## reaction happens
                    s[j].molecule_type = "NO2"
                    s[i].molecule_type = "O2"
                    set_new_velocites(s[i],s[j])
                    amount_of_reactions += 1
                    break
        elif s[j].molecule_type == "O3":
            for i in range(j+1, len(s)):
                if s[i].molecule_type == "NO" and reaction(s[j].body.velocity, s[i].body.velocity,
                            s[j].body.mass, s[i].body.mass):
                    ## reaction happens
                    s[j].molecule_type = "O2"
                    s[i].molecule_type = "NO2"
                    set_new_velocites(s[i],s[j])
                    amount_of_reactions += 1
                    break
    return True


def main():
    global amount_of_reactions

    ## construct space
    space = pymunk.Space()
    space.gravity = (0.0,0.0)
    ch = space.add_collision_handler(1, 1)
    ch.pre_solve = collision_handler
    addLines(space)

    ## list storing the amount of reactions for each round
    aor_list = np.zeros(MAX_ROUNDS)
    for r in range(1,MAX_ROUNDS):
        ## add new particles
        particles = addParticles(space, amount_list=[AMOUNT_OF_NO, AMOUNT_OF_O3], 
            molecule_type_list=["NO","O3"], existing_particles=[])
        ## display particles if GRAPHICS==True
        if GRAPHICS:
            pygame.init()
            screen = pygame.display.set_mode((SPACE_SIZE, SPACE_SIZE))
            pygame.display.set_caption("particles in gas")
            clock = pygame.time.Clock()

            draw_options = pymunk.pygame_util.DrawOptions(screen)

            for time in range(0, MAX_TIME):
                for event in pygame.event.get():
                    if event.type == QUIT:
                        sys.exit(0)
                    elif event.type == KEYDOWN and event.key == K_ESCAPE:
                        sys.exit(0)

                screen.fill((255,255,255))

                space.debug_draw(draw_options)

                space.step(DT)

                pygame.display.flip()
                clock.tick(1./DT)
        else:
            for time in range(0, MAX_TIME):
                space.step(DT)

        ## store and print amount of reactions in current round
        aor_list[r] = amount_of_reactions
        print("amount of reactions: {}".format(amount_of_reactions))
        
        ## delete all particles and reset amount of reactions
        for i in range (0,len(particles)):
            space.remove(list(particles[i].shapes)[0].body, list(particles[i].shapes)[0])
        particles=[]
        amount_of_reactions=0

    print("mean of first half:", np.mean(aor_list[range(1, int(MAX_ROUNDS/2))]))
    print("mean of second half:", np.mean(aor_list[range(int(MAX_ROUNDS/2), MAX_ROUNDS-1)]))

        

if __name__ == '__main__':
    main()