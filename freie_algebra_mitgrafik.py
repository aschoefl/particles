import sys, random
import pygame
from pygame.locals import *
import pymunk
import pymunk.pygame_util
import numpy as np 

random.seed = 42

## collision_type: 0-lines, 1-particle
## GLOBALS
SPACE_SIZE = 600
ELASTICITY = 0.999 # should not be >= 1 according to homepage
LINES_RADIUS = 4
AMOUNT_OF_NO = 300
AMOUNT_OF_O3 = 200
RADIUS_OF_PARTICLES = 3
MASS_OF_PARTICLES = 1
IMPULS_MEAN = 150
IMPULS_STD = 10
REACTION_VELOCITY = 100
MAX_TIME = 1000

amount_of_no2 = 0
amount_of_o2 = 0


def addParticles(space, amount, existing_particles = [], molecule_type=None):
    """Add a ball to the given space at a random position"""
    mass = MASS_OF_PARTICLES
    radius = RADIUS_OF_PARTICLES
    particles = existing_particles
    for j in range(0, amount):
        inertia = pymunk.moment_for_circle(mass, 0, radius, (0,0))
        body = pymunk.Body(mass, inertia)
        body.apply_impulse_at_local_point((random.choice([-1,1])*random.gauss(IMPULS_MEAN, IMPULS_STD), 
                                        random.choice([-1,1])*random.gauss(IMPULS_MEAN, IMPULS_STD)))
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

def collision_handler(arbiter, space, data):
    global amount_of_no2
    global amount_of_o2
    s = arbiter.shapes
    for j in range(0, len(s)):
        if s[j].molecule_type in ["NO2, O2"]:
            pass
        elif s[j].molecule_type == "NO":
            for i in range(j+1, len(s)):
                if s[i].molecule_type == "O3" and np.linalg.norm(s[j].body.velocity-s[i].body.velocity)>REACTION_VELOCITY:
                    s[j].molecule_type = "NO2"
                    s[i].molecule_type = "O2"
                    amount_of_no2+=1
                    amount_of_o2+=1
                    break
        elif s[j].molecule_type == "O3":
            for i in range(j+1, len(s)):
                if s[i].molecule_type == "NO" and np.linalg.norm(s[j].body.velocity-s[i].body.velocity)>REACTION_VELOCITY:
                    s[j].molecule_type = "O2"
                    s[i].molecule_type = "NO2"
                    amount_of_no2+=1
                    amount_of_o2+=1
                    break
                
    pass
    return True


def main():
    global amount_of_no2
    global amount_of_o2
    pygame.init()
    screen = pygame.display.set_mode((SPACE_SIZE, SPACE_SIZE))
    pygame.display.set_caption("Ich will nicht Algebra lernen")
    clock = pygame.time.Clock()

    space = pymunk.Space()
    space.gravity = (0.0,0.0)
    ch = space.add_collision_handler(1, 1)
    ch.pre_solve = collision_handler

    addLines(space)
    particles = addParticles(space, AMOUNT_OF_NO, molecule_type="NO")
    particles = addParticles(space, AMOUNT_OF_O3, existing_particles=particles, molecule_type="O3")
    draw_options = pymunk.pygame_util.DrawOptions(screen)

    for time in range(0, MAX_TIME):
        for event in pygame.event.get():
            if event.type == QUIT:
                sys.exit(0)
            elif event.type == KEYDOWN and event.key == K_ESCAPE:
                sys.exit(0)

        screen.fill((255,255,255))

        space.debug_draw(draw_options)

        space.step(1/50.0)

        pygame.display.flip()
        clock.tick(50)
    print("amount of no2: {} amount of o2: {}".format(amount_of_no2, amount_of_o2))

if __name__ == '__main__':
    main()