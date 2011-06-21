
*** Simulation.java ***

Represents the top most level of JAT simulation.  Create your run methods in here
using the SimModel.java class and it's methods.  Currently, there are several examples
provided.  The sequence is essentially this:

1 - Initialize jat.sim.SimModel
2 - Create your filenames and other relevant variables
3 - Create/initialize your spacecraft either from file or manually using jat.spacecraft.*
3.5 - this includes any controllers or state estimation (see jat.spacecraft.SpacecraftModel)
4 - Create/initialize your force models either from jat.forces.* or create your own
5 - Initialize the integration parameters: stepsize, thinning, t0, tf, MJD_UTC_Start
6 - Use the SimModel.runloop() method to propagate 
    OR manually propagate using SimModel.step(dt)

See jat.sim.SimModel for more details on simulation algorithms

To use Matlab functions, use the jat.matlabinterface.MatlabFunc class which allows Java to
call any Matlab function so long as JAT is run from the current Matlab session.

To use a Matlab controller, see jat.spacecraft.MatlabControlLaw

To use JAT from Matlab, use the "javaaddpath" command to add your JAT home directory and
the external JAR directories to Matlab.  Then, simply use the following syntax:

sim = jat.sim.SimModel()

force = jat.forces.ForceModel();

sim.addForce(force);

etc