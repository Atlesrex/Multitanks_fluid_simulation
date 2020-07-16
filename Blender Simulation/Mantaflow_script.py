
######################################################################
## LIBRARIES
######################################################################
from manta import *
import os.path, shutil, math, sys, gc, multiprocessing, platform, time

withMPBake = False # Bake files asynchronously
withMPSave = True # Save files asynchronously
isWindows = platform.system() != 'Darwin' and platform.system() != 'Linux'
# TODO (sebbas): Use this to simulate Windows multiprocessing (has default mode spawn)
#try:
#    multiprocessing.set_start_method('spawn')
#except:
#    pass

bpy = sys.modules.get('bpy')
if bpy is not None:
    sys.executable = bpy.app.binary_path_python

######################################################################
## VARIABLES
######################################################################

mantaMsg('Fluid variables')
dim_s7     = 3
res_s7     = 512
gravity_s7 = vec3(0, 0, -1)
gs_s7      = vec3(166, 54, 512)
maxVel_s7  = 0

doOpen_s7          = False
boundConditions_s7 = ''
boundaryWidth_s7   = 1

using_smoke_s7        = False
using_liquid_s7       = True
using_noise_s7        = False
using_adaptTime_s7    = True
using_obstacle_s7     = True
using_guiding_s7      = False
using_fractions_s7    = False
using_invel_s7        = False
using_outflow_s7      = False
using_sndparts_s7     = True
using_speedvectors_s7 = True

# Fluid time params
timeTotal_s7    = 0
timePerFrame_s7 = 0
frameLength_s7  = 0.0625
dt0_s7          = 0.0625
cflCond_s7      = 2
timestepsMin_s7 = 4
timestepsMax_s7 = 8

# Fluid diffusion / viscosity
domainSize_s7 = 1.7 # longest domain side in meters
viscosity_s7 = 1e-06 / (domainSize_s7*domainSize_s7) # kinematic viscosity in m^2/s

# Factor to convert blender velocities to manta velocities
toMantaUnitsFac_s7 = (1.0 / (1.0 / res_s7))
 # = dt/dx * 1/dt 
mantaMsg('Liquid variables')
narrowBandWidth_s7         = 3
combineBandWidth_s7        = narrowBandWidth_s7 - 1
adjustedNarrowBandWidth_s7 = 3 # only used in adjustNumber to control band width
particleNumber_s7   = 2
minParticles_s7     = 16
maxParticles_s7     = 32
radiusFactor_s7     = 1
using_mesh_s7       = True
using_final_mesh_s7 = True
using_fractions_s7  = False
fracThreshold_s7    = 0.05
flipRatio_s7        = 0.97
concaveUpper_s7     = 3.5
concaveLower_s7     = 0.4
meshRadiusFactor_s7 = 1
smoothenPos_s7      = 2
smoothenNeg_s7      = 2
randomness_s7       = 0.1
surfaceTension_s7   = 0

mantaMsg('Fluid variables mesh')
upres_sm7  = 2
gs_sm7     = vec3(upres_sm7*gs_s7.x, upres_sm7*gs_s7.y, upres_sm7*gs_s7.z)

mantaMsg('Fluid variables particles')
upres_sp7  = 2
gs_sp7     = vec3(upres_sp7*gs_s7.x, upres_sp7*gs_s7.y, upres_sp7*gs_s7.z)

tauMin_wc_sp7 = 2
tauMax_wc_sp7 = 8
tauMin_ta_sp7 = 5
tauMax_ta_sp7 = 20
tauMin_k_sp7 = 1
tauMax_k_sp7 = 5
k_wc_sp7 = 200
k_ta_sp7 = 60
k_b_sp7 = 0.5
k_d_sp7 = 0.6
lMin_sp7 = 12
lMax_sp7 = 25
c_s_sp7 = 0.4   # classification constant for snd parts
c_b_sp7 = 0.77  # classification constant for snd parts
pot_radius_sp7 = 3
update_radius_sp7 = 3
scaleFromManta_sp7 = 1.7 / float(res_s7) # resize factor for snd parts

######################################################################
## SOLVERS
######################################################################

mantaMsg('Solver base')
s7 = Solver(name='solver_base7', gridSize=gs_s7, dim=dim_s7)

mantaMsg('Solver mesh')
sm7 = Solver(name='solver_mesh7', gridSize=gs_sm7)

mantaMsg('Solver particles')
sp7 = Solver(name='solver_particles7', gridSize=gs_sp7)

######################################################################
## GRIDS
######################################################################

mantaMsg('Fluid alloc data')
flags_s7       = s7.create(FlagGrid)
vel_s7         = s7.create(MACGrid)
velC_s7        = s7.create(MACGrid)
x_vel_s7       = s7.create(RealGrid)
y_vel_s7       = s7.create(RealGrid)
z_vel_s7       = s7.create(RealGrid)
pressure_s7    = s7.create(RealGrid)
phiObs_s7      = s7.create(LevelsetGrid)
phiIn_s7       = s7.create(LevelsetGrid)
phiOut_s7      = s7.create(LevelsetGrid)
forces_s7      = s7.create(Vec3Grid)
x_force_s7     = s7.create(RealGrid)
y_force_s7     = s7.create(RealGrid)
z_force_s7     = s7.create(RealGrid)
obvel_s7       = None

# Keep track of important objects in dict to load them later on
fluid_data_dict_final_s7  = dict(vel=vel_s7)
fluid_data_dict_resume_s7 = dict(phiObs=phiObs_s7, phiIn=phiIn_s7, phiOut=phiOut_s7, flags=flags_s7)

mantaMsg('Liquid alloc')
phiParts_s7   = s7.create(LevelsetGrid)
phi_s7        = s7.create(LevelsetGrid)
phiTmp_s7     = s7.create(LevelsetGrid)
curvature_s7  = s7.create(RealGrid)
velOld_s7     = s7.create(MACGrid)
velParts_s7   = s7.create(MACGrid)
mapWeights_s7 = s7.create(MACGrid)
fractions_s7  = None # allocated dynamically

pp_s7         = s7.create(BasicParticleSystem)
pVel_pp7      = pp_s7.create(PdataVec3)

# Acceleration data for particle nbs
pindex_s7     = s7.create(ParticleIndexSystem)
gpi_s7        = s7.create(IntGrid)

# Keep track of important objects in dict to load them later on
liquid_data_dict_final_s7 = dict(pp=pp_s7, pVel=pVel_pp7)
liquid_data_dict_resume_s7 = dict(phiParts=phiParts_s7, phi=phi_s7, phiTmp=phiTmp_s7)

mantaMsg('Liquid alloc mesh')
phiParts_sm7 = sm7.create(LevelsetGrid)
phi_sm7      = sm7.create(LevelsetGrid)
pp_sm7       = sm7.create(BasicParticleSystem)
flags_sm7    = sm7.create(FlagGrid)
mesh_sm7     = sm7.create(Mesh)

if using_speedvectors_s7:
    mVel_mesh7 = mesh_sm7.create(MdataVec3)
    vel_sm7    = sm7.create(MACGrid)

# Acceleration data for particle nbs
pindex_sm7  = sm7.create(ParticleIndexSystem)
gpi_sm7     = sm7.create(IntGrid)

# Keep track of important objects in dict to load them later on
liquid_mesh_dict_s7 = dict(lMesh=mesh_sm7)

if using_speedvectors_s7:
    liquid_meshvel_dict_s7 = dict(lVelMesh=mVel_mesh7)

ppSnd_sp7         = sp7.create(BasicParticleSystem)
pVelSnd_pp7       = ppSnd_sp7.create(PdataVec3)
pForceSnd_pp7     = ppSnd_sp7.create(PdataVec3)
pLifeSnd_pp7      = ppSnd_sp7.create(PdataReal)
vel_sp7           = sp7.create(MACGrid)
flags_sp7         = sp7.create(FlagGrid)
phi_sp7           = sp7.create(LevelsetGrid)
phiIn_sp7         = sp7.create(LevelsetGrid)
phiObs_sp7        = sp7.create(LevelsetGrid)
phiObsIn_sp7      = sp7.create(LevelsetGrid)
normal_sp7        = sp7.create(VecGrid)
neighborRatio_sp7 = sp7.create(RealGrid)
trappedAir_sp7    = sp7.create(RealGrid)
waveCrest_sp7     = sp7.create(RealGrid)
kineticEnergy_sp7 = sp7.create(RealGrid)

# Keep track of important objects in dict to load them later on
liquid_particles_dict_final_s7  = dict(ppSnd=ppSnd_sp7, pVelSnd=pVelSnd_pp7, pLifeSnd=pLifeSnd_pp7)
liquid_particles_dict_resume_s7 = dict(trappedAir=trappedAir_sp7, waveCrest=waveCrest_sp7, kineticEnergy=kineticEnergy_sp7)

mantaMsg('Allocating obstacle data')
numObs_s7     = s7.create(RealGrid)
phiObsIn_s7   = s7.create(LevelsetGrid)
obvel_s7      = s7.create(MACGrid)
obvelC_s7     = s7.create(Vec3Grid)
x_obvel_s7    = s7.create(RealGrid)
y_obvel_s7    = s7.create(RealGrid)
z_obvel_s7    = s7.create(RealGrid)

if 'fluid_data_dict_resume_s7' in globals():
    fluid_data_dict_resume_s7.update(phiObsIn=phiObsIn_s7)

######################################################################
## DOMAIN INIT
######################################################################

# Prepare domain
phi_s7.initFromFlags(flags_s7)
phiIn_s7.initFromFlags(flags_s7)

######################################################################
## ADAPTIVE TIME
######################################################################

mantaMsg('Fluid adaptive time stepping')
s7.frameLength  = frameLength_s7
s7.timestepMin  = s7.frameLength / max(1, timestepsMax_s7)
s7.timestepMax  = s7.frameLength / max(1, timestepsMin_s7)
s7.cfl          = cflCond_s7
s7.timePerFrame = timePerFrame_s7
s7.timestep     = dt0_s7
s7.timeTotal    = timeTotal_s7
#mantaMsg('timestep: ' + str(s7.timestep) + ' // timPerFrame: ' + str(s7.timePerFrame) + ' // frameLength: ' + str(s7.frameLength) + ' // timeTotal: ' + str(s7.timeTotal) )

def fluid_adapt_time_step_7():
    mantaMsg('Fluid adapt time step')
    
    # time params are animatable
    s7.frameLength = frameLength_s7
    s7.cfl         = cflCond_s7
    
    # ensure that vel grid is full (remember: adaptive domain can reallocate solver)
    copyRealToVec3(sourceX=x_vel_s7, sourceY=y_vel_s7, sourceZ=z_vel_s7, target=vel_s7)
    maxVel_s7 = vel_s7.getMax() if vel_s7 else 0
    if using_adaptTime_s7:
        mantaMsg('Adapt timestep, maxvel: ' + str(maxVel_s7))
        s7.adaptTimestep(maxVel_s7)

######################################################################
## IMPORT
######################################################################

def fluid_file_import_s7(dict, path, framenr, file_format):
    try:
        framenr = fluid_cache_get_framenr_formatted_7(framenr)
        for name, object in dict.items():
            file = os.path.join(path, name + '_' + framenr + file_format)
            if os.path.isfile(file):
                object.load(file)
            else:
                mantaMsg('Could not load file ' + str(file))
    except:
        mantaMsg('exception found')
        #mantaMsg(str(e))
        pass # Just skip file load errors for now

def fluid_cache_get_framenr_formatted_7(framenr):
    return str(framenr).zfill(4) # framenr with leading zeroes

def fluid_load_data_7(path, framenr, file_format, resumable):
    mantaMsg('Fluid load data, frame ' + str(framenr))
    fluid_file_import_s7(dict=fluid_data_dict_final_s7, path=path, framenr=framenr, file_format=file_format)
    
    if resumable:
        fluid_file_import_s7(dict=fluid_data_dict_resume_s7, path=path, framenr=framenr, file_format=file_format)
        
        # When adaptive domain bake is resumed we need correct values in xyz vel grids
        copyVec3ToReal(source=vel_s7, targetX=x_vel_s7, targetY=y_vel_s7, targetZ=z_vel_s7)

def liquid_load_data_7(path, framenr, file_format, resumable):
    mantaMsg('Liquid load data')
    fluid_file_import_s7(dict=liquid_data_dict_final_s7, path=path, framenr=framenr, file_format=file_format)
    if resumable:
        fluid_file_import_s7(dict=liquid_data_dict_resume_s7, path=path, framenr=framenr, file_format=file_format)

def liquid_load_mesh_7(path, framenr, file_format):
    mantaMsg('Liquid load mesh')
    fluid_file_import_s7(dict=liquid_mesh_dict_s7, path=path, framenr=framenr, file_format=file_format)

def liquid_load_meshvel_7(path, framenr, file_format):
    mantaMsg('Liquid load meshvel')
    fluid_file_import_s7(dict=liquid_meshvel_dict_s7, path=path, framenr=framenr, file_format=file_format)

def liquid_load_particles_7(path, framenr, file_format, resumable):
    mantaMsg('Liquid load particles')
    fluid_file_import_s7(dict=liquid_particles_dict_final_s7, path=path, framenr=framenr, file_format=file_format)
    if resumable:
        fluid_file_import_s7(dict=liquid_particles_dict_resume_s7, path=path, framenr=framenr, file_format=file_format)

######################################################################
## PRE/POST STEPS
######################################################################

def fluid_pre_step_7():
    mantaMsg('Fluid pre step')
    
    phiObs_s7.setConst(9999)
    phiOut_s7.setConst(9999)
    
    # Main vel grid is copied in adapt time step function
    
    # translate obvels (world space) to grid space
    if using_obstacle_s7:
        # Average out velocities from multiple obstacle objects at one cell
        x_obvel_s7.safeDivide(numObs_s7)
        y_obvel_s7.safeDivide(numObs_s7)
        z_obvel_s7.safeDivide(numObs_s7)
        
        x_obvel_s7.multConst(toMantaUnitsFac_s7)
        y_obvel_s7.multConst(toMantaUnitsFac_s7)
        z_obvel_s7.multConst(toMantaUnitsFac_s7)
        
        copyRealToVec3(sourceX=x_obvel_s7, sourceY=y_obvel_s7, sourceZ=z_obvel_s7, target=obvelC_s7)
    
    # translate invels (world space) to grid space
    if using_invel_s7:
        x_invel_s7.multConst(toMantaUnitsFac_s7)
        y_invel_s7.multConst(toMantaUnitsFac_s7)
        z_invel_s7.multConst(toMantaUnitsFac_s7)
        copyRealToVec3(sourceX=x_invel_s7, sourceY=y_invel_s7, sourceZ=z_invel_s7, target=invelC_s7)
    
    if using_guiding_s7:
        weightGuide_s7.multConst(0)
        weightGuide_s7.addConst(alpha_sg7)
        interpolateMACGrid(source=guidevel_sg7, target=velT_s7)
        velT_s7.multConst(vec3(gamma_sg7))
    
    # translate external forces (world space) to grid space
    x_force_s7.multConst(toMantaUnitsFac_s7)
    y_force_s7.multConst(toMantaUnitsFac_s7)
    z_force_s7.multConst(toMantaUnitsFac_s7)
    copyRealToVec3(sourceX=x_force_s7, sourceY=y_force_s7, sourceZ=z_force_s7, target=forces_s7)
    
    # If obstacle has velocity, i.e. is a moving obstacle, switch to dynamic preconditioner
    if using_smoke_s7 and using_obstacle_s7 and obvelC_s7.getMax() > 0:
        mantaMsg('Using dynamic preconditioner')
        preconditioner_s7 = PcMGDynamic
    else:
        mantaMsg('Using static preconditioner')
        preconditioner_s7 = PcMGStatic

def fluid_post_step_7():
    mantaMsg('Fluid post step')
    forces_s7.clear()
    x_force_s7.clear()
    y_force_s7.clear()
    z_force_s7.clear()
    
    if using_guiding_s7:
        weightGuide_s7.clear()
    if using_invel_s7:
        x_invel_s7.clear()
        y_invel_s7.clear()
        z_invel_s7.clear()
        invel_s7.clear()
        invelC_s7.clear()
    
    # Copy vel grid to reals grids (which Blender internal will in turn use for vel access)
    copyVec3ToReal(source=vel_s7, targetX=x_vel_s7, targetY=y_vel_s7, targetZ=z_vel_s7)

######################################################################
## STEPS
######################################################################

def liquid_adaptive_step_7(framenr):
    mantaMsg('Manta step, frame ' + str(framenr))
    s7.frame = framenr
    
    fluid_pre_step_7()
    
    flags_s7.initDomain(boundaryWidth=1 if using_fractions_s7 else 0, phiWalls=phiObs_s7, outflow=boundConditions_s7)
    
    if using_obstacle_s7:
        mantaMsg('Initializing obstacle levelset')
        phiObsIn_s7.fillHoles(maxDepth=int(res_s7), boundaryWidth=2)
        extrapolateLsSimple(phi=phiObsIn_s7, distance=int(res_s7/2), inside=True)
        extrapolateLsSimple(phi=phiObsIn_s7, distance=3, inside=False)
        phiObs_s7.join(phiObsIn_s7)
        
        # Using boundaryWidth=2 to not search beginning from walls (just a performance optimization)
        # Additional sanity check: fill holes in phiObs which can result after joining with phiObsIn
        phiObs_s7.fillHoles(maxDepth=int(res_s7), boundaryWidth=2)
        extrapolateLsSimple(phi=phiObs_s7, distance=int(res_s7/2), inside=True)
        extrapolateLsSimple(phi=phiObs_s7, distance=3, inside=False)
    
    mantaMsg('Initializing fluid levelset')
    extrapolateLsSimple(phi=phiIn_s7, distance=int(res_s7/2), inside=True)
    extrapolateLsSimple(phi=phiIn_s7, distance=int(res_s7/2), inside=False)
    phi_s7.join(phiIn_s7)
    
    if using_obstacle_s7:
        phi_s7.subtract(o=phiObsIn_s7, flags=flags_s7, subtractType=FlagObstacle)
    
    if using_outflow_s7:
        phiOut_s7.join(phiOutIn_s7)
    
    if using_fractions_s7:
        updateFractions(flags=flags_s7, phiObs=phiObs_s7, fractions=fractions_s7, boundaryWidth=boundaryWidth_s7, fracThreshold=fracThreshold_s7)
    setObstacleFlags(flags=flags_s7, phiObs=phiObs_s7, phiOut=phiOut_s7, fractions=fractions_s7, phiIn=phiIn_s7)
    
    # add initial velocity: set invel as source grid to ensure const vels in inflow region, sampling makes use of this
    if using_invel_s7:
        extrapolateVec3Simple(vel=invelC_s7, phi=phiIn_s7, distance=int(res_s7/2), inside=True)
        resampleVec3ToMac(source=invelC_s7, target=invel_s7)
        pVel_pp7.setSource(invel_s7, isMAC=True)
    
    sampleLevelsetWithParticles(phi=phiIn_s7, flags=flags_s7, parts=pp_s7, discretization=particleNumber_s7, randomness=randomness_s7)
    flags_s7.updateFromLevelset(phi_s7)
    
    mantaMsg('Liquid step / s7.frame: ' + str(s7.frame))
    liquid_step_7()
    
    s7.step()
    
    fluid_post_step_7()

def liquid_step_7():
    mantaMsg('Liquid step')
    
    mantaMsg('Advecting particles')
    pp_s7.advectInGrid(flags=flags_s7, vel=vel_s7, integrationMode=IntRK4, deleteInObstacle=False, stopInObstacle=False)
    
    mantaMsg('Pushing particles out of obstacles')
    pushOutofObs(parts=pp_s7, flags=flags_s7, phiObs=phiObs_s7)
    
    mantaMsg('Advecting phi')
    advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=phi_s7, order=1) # first order is usually enough
    mantaMsg('Advecting velocity')
    advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=vel_s7, order=2)
    
    phiTmp_s7.copyFrom(phi_s7) # save original phi for later use in mesh creation
    
    # create level set of particles
    gridParticleIndex(parts=pp_s7, flags=flags_s7, indexSys=pindex_s7, index=gpi_s7)
    unionParticleLevelset(parts=pp_s7, indexSys=pindex_s7, flags=flags_s7, index=gpi_s7, phi=phiParts_s7, radiusFactor=radiusFactor_s7)
    
    # combine level set of particles with grid level set
    phi_s7.addConst(1.) # shrink slightly
    phi_s7.join(phiParts_s7)
    extrapolateLsSimple(phi=phi_s7, distance=narrowBandWidth_s7+2, inside=True)
    extrapolateLsSimple(phi=phi_s7, distance=3)
    phi_s7.setBoundNeumann(0) # make sure no particles are placed at outer boundary
    
    if doOpen_s7 or using_outflow_s7:
        resetOutflow(flags=flags_s7, phi=phi_s7, parts=pp_s7, index=gpi_s7, indexSys=pindex_s7)
    flags_s7.updateFromLevelset(phi_s7)
    
    # combine particles velocities with advected grid velocities
    mapPartsToMAC(vel=velParts_s7, flags=flags_s7, velOld=velOld_s7, parts=pp_s7, partVel=pVel_pp7, weight=mapWeights_s7)
    extrapolateMACFromWeight(vel=velParts_s7, distance=2, weight=mapWeights_s7)
    combineGridVel(vel=velParts_s7, weight=mapWeights_s7, combineVel=vel_s7, phi=phi_s7, narrowBand=combineBandWidth_s7, thresh=0)
    velOld_s7.copyFrom(vel_s7)
    
    # forces & pressure solve
    addGravity(flags=flags_s7, vel=vel_s7, gravity=gravity_s7)
    
    mantaMsg('Adding external forces')
    addForceField(flags=flags_s7, vel=vel_s7, force=forces_s7)
    
    if using_obstacle_s7:
        mantaMsg('Extrapolating object velocity')
        # ensure velocities inside of obs object, slightly add obvels outside of obs object
        extrapolateVec3Simple(vel=obvelC_s7, phi=phiObsIn_s7, distance=int(res_s7/2), inside=True)
        extrapolateVec3Simple(vel=obvelC_s7, phi=phiObsIn_s7, distance=3, inside=False)
        resampleVec3ToMac(source=obvelC_s7, target=obvel_s7)
    
    extrapolateMACSimple(flags=flags_s7, vel=vel_s7, distance=2, intoObs=True if using_fractions_s7 else False)
    
    # vel diffusion / viscosity!
    if viscosity_s7 > 0.:
        mantaMsg('Viscosity')
        # diffusion param for solve = const * dt / dx^2
        alphaV = viscosity_s7 * s7.timestep * float(res_s7*res_s7)
        setWallBcs(flags=flags_s7, vel=vel_s7, obvel=None if using_fractions_s7 else obvel_s7, phiObs=phiObs_s7, fractions=fractions_s7)
        cgSolveDiffusion(flags_s7, vel_s7, alphaV)
    
    setWallBcs(flags=flags_s7, vel=vel_s7, obvel=None if using_fractions_s7 else obvel_s7, phiObs=phiObs_s7, fractions=fractions_s7)
    
    mantaMsg('Calculating curvature')
    getLaplacian(laplacian=curvature_s7, grid=phi_s7)
    
    if using_guiding_s7:
        mantaMsg('Guiding and pressure')
        PD_fluid_guiding(vel=vel_s7, velT=velT_s7, flags=flags_s7, phi=phi_s7, curv=curvature_s7, surfTens=surfaceTension_s7, fractions=fractions_s7, weight=weightGuide_s7, blurRadius=beta_sg7, pressure=pressure_s7, tau=tau_sg7, sigma=sigma_sg7, theta=theta_sg7, zeroPressureFixing=not doOpen_s7)
    else:
        mantaMsg('Pressure')
        solvePressure(flags=flags_s7, vel=vel_s7, pressure=pressure_s7, phi=phi_s7, curv=curvature_s7, surfTens=surfaceTension_s7, fractions=fractions_s7, obvel=obvel_s7 if using_fractions_s7 else None)
    
    extrapolateMACSimple(flags=flags_s7, vel=vel_s7, distance=4, intoObs=True if using_fractions_s7 else False)
    setWallBcs(flags=flags_s7, vel=vel_s7, obvel=None if using_fractions_s7 else obvel_s7, phiObs=phiObs_s7, fractions=fractions_s7)
    
    if not using_fractions_s7:
        extrapolateMACSimple(flags=flags_s7, vel=vel_s7, distance=(int(maxVel_s7*1.25)))
    
    # set source grids for resampling, used in adjustNumber!
    pVel_pp7.setSource(vel_s7, isMAC=True)
    adjustNumber(parts=pp_s7, vel=vel_s7, flags=flags_s7, minParticles=minParticles_s7, maxParticles=maxParticles_s7, phi=phi_s7, exclude=phiObs_s7, radiusFactor=radiusFactor_s7, narrowBand=adjustedNarrowBandWidth_s7)
    flipVelocityUpdate(vel=vel_s7, velOld=velOld_s7, flags=flags_s7, parts=pp_s7, partVel=pVel_pp7, flipRatio=flipRatio_s7)

def liquid_step_mesh_7():
    mantaMsg('Liquid step mesh')
    
    interpolateGrid(target=phi_sm7, source=phiTmp_s7)
    
    # create surface
    pp_sm7.readParticles(pp_s7)
    gridParticleIndex(parts=pp_sm7, flags=flags_sm7, indexSys=pindex_sm7, index=gpi_sm7)
    
    if using_final_mesh_s7:
        mantaMsg('Liquid using improved particle levelset')
        improvedParticleLevelset(pp_sm7, pindex_sm7, flags_sm7, gpi_sm7, phiParts_sm7, meshRadiusFactor_s7, smoothenPos_s7, smoothenNeg_s7, concaveLower_s7, concaveUpper_s7)
    else:
        mantaMsg('Liquid using union particle levelset')
        unionParticleLevelset(pp_sm7, pindex_sm7, flags_sm7, gpi_sm7, phiParts_sm7, meshRadiusFactor_s7)
    
    phi_sm7.addConst(1.) # shrink slightly
    phi_sm7.join(phiParts_sm7)
    extrapolateLsSimple(phi=phi_sm7, distance=narrowBandWidth_s7+2, inside=True)
    extrapolateLsSimple(phi=phi_sm7, distance=3)
    phi_sm7.setBoundNeumann(0) # make sure no particles are placed at outer boundary
    
    # Vert vel vector needs to pull data from vel grid with correct dim
    if using_speedvectors_s7:
        interpolateMACGrid(target=vel_sm7, source=vel_s7)
        mVel_mesh7.setSource(vel_sm7, isMAC=True)
    
    phi_sm7.setBound(0.5,int(((upres_sm7)*2)-2) )
    phi_sm7.createMesh(mesh_sm7)

def liquid_step_particles_7():
    mantaMsg('Secondary particles step')
    
    # no upres: just use the loaded grids
    if upres_sp7 <= 1:
        flags_sp7.copyFrom(flags_s7)
        vel_sp7.copyFrom(vel_s7)
        phiObs_sp7.copyFrom(phiObs_s7)
        phi_sp7.copyFrom(phi_s7)
    
    # with upres: recreate grids
    else:
        # create highres grids by interpolation
        interpolateMACGrid(target=vel_sp7, source=vel_s7)
        interpolateGrid(target=phi_sp7, source=phi_s7)
        flags_sp7.initDomain(boundaryWidth=boundaryWidth_s7, phiWalls=phiObs_sp7, outflow=boundConditions_s7)
        flags_sp7.updateFromLevelset(levelset=phi_sp7)
    
    # actual secondary simulation
    #extrapolateLsSimple(phi=phi_sp7, distance=radius+1, inside=True)
    flipComputeSecondaryParticlePotentials(potTA=trappedAir_sp7, potWC=waveCrest_sp7, potKE=kineticEnergy_sp7, neighborRatio=neighborRatio_sp7, flags=flags_sp7, v=vel_sp7, normal=normal_sp7, phi=phi_sp7, radius=pot_radius_sp7, tauMinTA=tauMin_ta_sp7, tauMaxTA=tauMax_ta_sp7, tauMinWC=tauMin_wc_sp7, tauMaxWC=tauMax_wc_sp7, tauMinKE=tauMin_k_sp7, tauMaxKE=tauMax_k_sp7, scaleFromManta=scaleFromManta_sp7)
    flipSampleSecondaryParticles(mode='single', flags=flags_sp7, v=vel_sp7, pts_sec=ppSnd_sp7, v_sec=pVelSnd_pp7, l_sec=pLifeSnd_pp7, lMin=lMin_sp7, lMax=lMax_sp7, potTA=trappedAir_sp7, potWC=waveCrest_sp7, potKE=kineticEnergy_sp7, neighborRatio=neighborRatio_sp7, c_s=c_s_sp7, c_b=c_b_sp7, k_ta=k_ta_sp7, k_wc=k_wc_sp7, dt=s7.frameLength)
    flipUpdateSecondaryParticles(mode='linear', pts_sec=ppSnd_sp7, v_sec=pVelSnd_pp7, l_sec=pLifeSnd_pp7, f_sec=pForceSnd_pp7, flags=flags_sp7, v=vel_sp7, neighborRatio=neighborRatio_sp7, radius=update_radius_sp7, gravity=gravity_s7, k_b=k_b_sp7, k_d=k_d_sp7, c_s=c_s_sp7, c_b=c_b_sp7, dt=s7.frameLength)
    if 1:
        pushOutofObs(parts = ppSnd_sp7, flags = flags_sp7, phiObs = phiObs_sp7, shift = 1.0)
    flipDeleteParticlesInObstacle(pts=ppSnd_sp7, flags=flags_sp7)
    #debugGridInfo(flags = flags_sp7, grid = trappedAir_sp7, name = 'Trapped Air')
    #debugGridInfo(flags = flags_sp7, grid = waveCrest_sp7, name = 'Wave Crest')
    #debugGridInfo(flags = flags_sp7, grid = kineticEnergy_sp7, name = 'Kinetic Energy')

######################################################################
## MAIN
######################################################################

# Helper function to call cache load functions
def load(frame, cache_resumable):
    fluid_load_data_7(os.path.join(cache_dir, 'data'), frame, file_format_data, cache_resumable)
    liquid_load_data_7(os.path.join(cache_dir, 'data'), frame, file_format_data, cache_resumable)
    if using_sndparts_s7:
        liquid_load_particles_7(os.path.join(cache_dir, 'particles'), frame, file_format_particles, cache_resumable)
    if using_mesh_s7:
        liquid_load_mesh_7(os.path.join(cache_dir, 'mesh'), frame, file_format_mesh)
    if using_guiding_s7:
        fluid_load_guiding_7(os.path.join(cache_dir, 'guiding'), frame, file_format_data)

# Helper function to call step functions
def step(frame):
    liquid_adaptive_step_7(frame)
    if using_mesh_s7:
        liquid_step_mesh_7()
    if using_sndparts_s7:
        liquid_step_particles_7()

gui = None
if (GUI):
    gui=Gui()
    gui.show()
    gui.pause()

cache_resumable       = True
cache_dir             = 'C:\Users\Aleksei\Google Drive\__current_data\cache_fluid'
file_format_data      = '.uni'
file_format_noise     = '.uni'
file_format_particles = '.uni'
file_format_mesh      = '.bobj.gz'

# Start and stop for simulation
current_frame  = 1
end_frame      = 240

# How many frame to load from cache
from_cache_cnt = 100

loop_cnt = 0
while current_frame <= end_frame:
    
    # Load already simulated data from cache:
    if loop_cnt < from_cache_cnt:
        load(current_frame, cache_resumable)
    
    # Otherwise simulate new data
    else:
        while(s7.frame <= current_frame):
            if using_adaptTime_s7:
                fluid_adapt_time_step_7()
            step(current_frame)
    
    current_frame += 1
    loop_cnt += 1
    
    if gui:
        gui.pause()
