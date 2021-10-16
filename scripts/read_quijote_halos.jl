using PyCall

homedir = @__DIR__

py""" 
import imp
"""

readsnap = py"imp.load_source"("readsnap",string(homedir,"/readsnap.py"))
readgadget = py"imp.load_source"("readgadget", string(homedir,"/readgadget.py"))
readfof = py"imp.load_source"("readfof",string(homedir,"/readfof.py"))

"""
Return information about halos.

Input:
- dirname: string, name of the directory
- redshift_int: int, 4 -> 0.0, ... , 0 -> 3.0, (see QUIJOTE documentation)
Output:
- halos: has fields like halos.GroupPos, halos.GroupVel, halos.GroupMass.
GroupVel needs to be multiplied by 1 + z.
"""
function read_halos(dirname, redshift_int)
    halos = readfof.FoF_catalog(dirname, redshift_int, long_ids=false, swap=false, SFR=false, read_IDs=false)
end
