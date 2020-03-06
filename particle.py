#!python3
import sys, random
import pymunk
import numpy as np 
from scipy import constants

import pygame
from pygame.locals import *
import pymunk.pygame_util
 

random.seed = 42

## collision_type: 0-lines, 1-particle

## GLOBALS
GRAPHICS = True
T = 273 # Kelvin
k = 1.38064852 # eigentlich 10^-23
SPACE_SIZE = 1000
ELASTICITY = 0.999 # should not be >= 1 according to docu of pymunk
LINES_RADIUS = 4
AMOUNT_OF_NO = 300
AMOUNT_OF_O3 = 300
RADIUS_OF_PARTICLES = 1 # eigentlich 10^-10 m
MASS_OF_NO = (2.3+2.67)  # eigentlich 10^-26 kg
MASS_OF_O3 = 2.67*3     # eigentlich 10^-26 kg
VELOCITY_MEAN_NO = (4*k*T/MASS_OF_NO)**(1/2)*10**(3/2)
VELOCITY_STD_NO = (k*T/MASS_OF_NO*(2-(2**(3/2)-4))*constants.pi)**(1/2)*10**(3/2)  #(k*T/MASS_OF_NO*(4*constants.pi-2))**(1/2)
VELOCITY_MEAN_O3 = (4*k*T/MASS_OF_O3)**(1/2)*10**(3/2)
VELOCITY_STD_O3 = (k*T/MASS_OF_NO*(2-(2**(3/2)-4))*constants.pi)**(1/2)*10**(3/2)  #(k*T/MASS_OF_O3*(4*constants.pi-2))**(1/2)
REACTION_ENERGY = 1.93*10**6 # eigentlich 10^-20 J
MAX_TIME = 500
MAX_ROUNDS = 300
DT = 0.01 # timestep was 0.001

amount_of_reactions = 0


def addParticles(space, amount_list, molecule_type_list, existing_particles = []):
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
            p_abs = random.gauss(vel_mean, vel_std)*mass
            x_dir = random.random()*random.choice([-1,1])
            y_dir = (1-x_dir**2)**(1/2)*random.choice([-1,1])
            body.apply_impulse_at_local_point((x_dir*p_abs, y_dir*p_abs))
            x = random.uniform(radius/2+LINES_RADIUS/2+2,SPACE_SIZE-radius/2-LINES_RADIUS/2-2)
            y = random.uniform(radius/2+LINES_RADIUS/2+2,SPACE_SIZE-radius/2-LINES_RADIUS/2-2)
            body.position = x, y
            shape = pymunk.Circle(body, radius, (0,0))
            shape.elasticity = ELASTICITY
            shape.friction = 0
            shape.collision_type = 1
            shape.molecule_type = molecule_type
            space.add(body, shape)
            try:
                particles.append(body)
            except:
                "Wrong type for existing_particles in addParticles"
    return(particles)


def addLines(space):
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
        l.collision_type = 0
        
    space.add(lines[0], lines[1], lines[2], lines[3], body)

def reaction(v1, v2, m1, m2):
    try:
        if REACTION_ENERGY*2*(m1+m2)/((v1-v2)**2*m1*m2) >= 0:
            return True
        else:
            return False
    except:
        return False

def collision_handler(arbiter, space, data):
    global amount_of_reactions
    s = arbiter.shapes
    n = np.array(arbiter.contact_point_set.normal)
    for j in range(0, len(s)):
        if s[j].molecule_type in ["NO2, O2"]:
            pass
        elif s[j].molecule_type == "NO":
            for i in range(j+1, len(s)):
                if s[i].molecule_type == "O3" and reaction(np.dot(np.array(s[j].body.velocity),n), np.dot(np.array(s[i].body.velocity),n),
                            s[j].body.mass, s[i].body.mass):
                    s[j].molecule_type = "NO2"
                    s[i].molecule_type = "O2"
                    amount_of_reactions += 1
                    #print("reaction")
                    break
        elif s[j].molecule_type == "O3":
            for i in range(j+1, len(s)):
                if s[i].molecule_type == "NO" and reaction(np.dot(np.array(s[j].body.velocity),n), np.dot(np.array(s[i].body.velocity),n),
                            s[j].body.mass, s[i].body.mass):
                    s[j].molecule_type = "O2"
                    s[i].molecule_type = "NO2"
                    amount_of_reactions += 1
                    #print("reaction")
                    break
                
    pass
    return True


def main():
    global amount_of_reactions

    space = pymunk.Space()
    space.gravity = (0.0,0.0)
    ch = space.add_collision_handler(1, 1)
    ch.pre_solve = collision_handler

    addLines(space)
    aor_list = np.zeros(MAX_ROUNDS)
    for r in range(1,MAX_ROUNDS):
        particles = addParticles(space, amount_list=[AMOUNT_OF_NO, AMOUNT_OF_O3], molecule_type_list=["NO","O3"], existing_particles=[])
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

                space.step(DT) # was 1./50

                pygame.display.flip()
                clock.tick(1./DT)
        else:
            for time in range(0, MAX_TIME):
                space.step(DT)
        aor_list[r] = amount_of_reactions
        print("amount of reactions: {}".format(amount_of_reactions))
        
        for i in range (0,len(particles)):
            space.remove(list(particles[i].shapes)[0].body, list(particles[i].shapes)[0])
        particles=[]
        
        amount_of_reactions=0

    print("mean of first half:", np.mean(aor_list[range(1, int(MAX_ROUNDS/2))]))
    print("mean of second half:", np.mean(aor_list[range(int(MAX_ROUNDS/2), MAX_ROUNDS-1)]))

        

if __name__ == '__main__':
    main()