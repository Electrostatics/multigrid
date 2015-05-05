#Types of Multigrid
There are currently four mg calculations in APBS: mg-auto, mg-dummy, mg-manual, and mg-para.

##mg-auto
This multigrid calculation automatically sets up and performs a string of single-point PBE calculations to "focus" on a region of interest (binding site, etc.) in a system.  It is basically an automated version of mg-manual designed for easier use.

##mg-dummy
This type of calculation allows users to write out dielectric, ion-accessibility, and charge distribution, and other types of maps that depend solely on biomolecular geometry. Since these maps depend only on geometry, they can be written out without actually solving the PB equation. The syntax for this command is identical to mg-manual.

This one is used to setup the problem, and can save time if running different calculations on the same biomolecule.

##mg-manual
This is a standard single-point multigrid PBE calculation without focusing or additional refinement. The mg-manual calculation offers the most control of parameters to the user. Several of these calculations can be strung together to perform focusing calculations by judicious choice of the bcfl flag; however, the setup of the focusing is not automated as it is in mg-auto and mg-para calculations and therefore this command should primarily be used by more experienced users.

##mg-para
This calculation closely resembles mg-auto in syntax. However, it is designed to perform electrostatics calculations on systems in a parallel focusing fashion.

#Parameters to Multigrid
Name|Parameters|Required?|auto|dummy|manual|para
---:|-----|:-------:|:--:|:---:|:----:|:--:
async|`integer`|||||✅
bcfl|`zero`, `sdh`, `mdh`, `focus`, `map`|✅|✅|✅|✅|✅
calcenergy|`no`, `total`, `comps`||✅|✅|✅|✅
calcforce|`no`, `total`, `comps`||✅|✅|✅|✅
chgm |`sp10`, `sp12`, `sp14`|✅|✅|✅|✅|✅
dime |`integer integer integer`|✅|✅|||✅
etol |`float`||✅|✅|✅|✅
cgcent|`mol {molecule_id}`, `float float float`|✅|✅|✅|✅|✅
cglen|`float float float`|✅|✅|✅|✅|✅
fgcent|`mol {molecule_id}`, `float float float`|✅|✅|✅|✅|✅
fglen|`float float float`|✅|✅|✅|✅|✅
gcent|`mol {molecule_id}`, `float float float`|✅||✅|✅|✅
glen|`float float float`|✅||✅|✅|✅
grid|`float float float`|✅||✅|✅|✅
ion|`charge {float} conc {float} radius {float}`||✅|✅|✅|✅
lpbe||🌟(npbe)|✅|✅|✅|✅
mol|`molecule_id`|✅|✅|✅|✅|✅
nlev|`integer`|✅||✅|✅|✅
npbe||🌟(lpbe)|✅|✅|✅|✅
pdie|`float`|✅|✅|✅|✅|✅
pdime|`integer integer integer`|✅||||✅
sdens|`float`|✅|✅|✅|✅|✅
sdie|`float`|✅|✅|✅|✅|✅
smpbe|`vol {float} size {float}`||✅|✅|✅|✅
srad|`float`|🌟(swin)|✅|✅|✅|✅
srfm|`mol`, `smol`, `spl2`, `spl4`|✅|✅|✅|✅|✅
swin|`float`|🌟(srad)|✅|✅|✅|✅
temp|`float`|✅|✅|✅|✅|✅
useaqua|||✅|✅|✅|✅
usemap|`type [diel, kappa, charge, pot] {file_id}`|✅|✅|✅|✅|✅
usemesh|`{file_id}`|✅|✅|✅|✅|✅
write|`type format stem` where *type* is one of `charge`, `pot`, `atompot`, `smol`, `sspl`, `vdw`, `ivdw`, `lap`, `edens`, `ndens`, `qdens`, `dielx`, `diely`, `dielz`, or `kappa`, *format* is one of `dx`, `gz`, `uhbd`, or `flat`, and *stem* is an output path|✅|✅|✅|✅|✅
writemat|`poisson {output path}`|✅|✅|✅|✅|✅

#Datastructures
##Atom Representation
```C
struct sVatom {

    double position[3];  /**< Atomic position */
    double radius;  /**< Atomic radius   */
    double charge;  /**< Atomic charge   */
    double partID;  /**< Partition value for assigning atoms to particular
                     * processors and/or partitions   */
    double epsilon; /**< Epsilon value for WCA calculations */

    int id;  /**< Atomic ID; this should be a unique non-negative integer
              * assigned based on the index of the atom in a Valist atom
              * array */

    char resName[VMAX_RECLEN]; /**< Residue name from PDB/PQR file */
    char atomName[VMAX_RECLEN]; /**< Atom name from PDB/PDR file */

#if defined(WITH_TINKER)

    double dipole[3];          /**< Permanent dipole */
    double quadrupole[9];      /**< Permanent quadrupole */
    double inducedDipole[3];   /**< Induced dipole */
    double nlInducedDipole[3];  /**< Non-local induced dipole */

#endif /* if defined(WITH_TINKER) */
};
```

##Input File Parsing
```C
struct sNOsh {

    NOsh_calc *calc[NOSH_MAXCALC];  /**< The array of calculation objects
        corresponding to actual calculations performed by the code.  Compare to
        sNOsh::elec */
    int ncalc;  /**< The number of calculations in the calc array */

    NOsh_calc *elec[NOSH_MAXCALC];  /**< The array of calculation objects
        corresponding to ELEC statements read in the input file.  Compare to
        sNOsh::calc */
    int nelec;  /**< The number of elec statements in the input file and in the
        elec array */

    NOsh_calc *apol[NOSH_MAXCALC];  /**< The array of calculation objects
        corresponding to APOLAR statements read in the input file.  Compare to
        sNOsh::calc */
    int napol;  /**< The number of apolar statements in the input file and in the
        apolar array */

    int ispara;  /**< 1 => is a parallel calculation, 0 => is not */
    int proc_rank;  /**< Processor rank in parallel calculation */
    int proc_size;  /**< Number of processors in parallel calculation */
    int bogus;  /**< A flag which tells routines using NOsh that this particular
        NOsh is broken -- useful for parallel focusing calculations where the
        user gave us too many processors (1 => ignore this NOsh; 0 => this NOsh
                                          is OK) */
    int elec2calc[NOSH_MAXCALC];  /**< A mapping between ELEC statements which
        appear in the input file and calc objects stored above.  Since we allow
        both normal and focused  multigrid, there isn't a 1-to-1 correspondence
        between ELEC statements and actual calcualtions.  This can really
        confuse operations which work on specific calculations further down the
        road (like PRINT).  Therefore this array is the initial point of entry
        for any calculation-specific operation.  It points to a specific entry
        in the calc array. */
    int apol2calc[NOSH_MAXCALC];  /**< (see elec2calc) */

    int nmol;  /**< Number of molecules */
    char molpath[NOSH_MAXMOL][VMAX_ARGLEN];   /**< Paths to mol files */
    NOsh_MolFormat molfmt[NOSH_MAXMOL];  /**< Mol files formats */
    Valist *alist[NOSH_MAXMOL];  /**<  Molecules for calculation (can be used in
        setting mesh centers */
    int gotparm;  /**< Either have (1) or don't have (0) parm */
    char parmpath[VMAX_ARGLEN];   /**< Paths to parm file */
    NOsh_ParmFormat parmfmt;  /**< Parm file format */
    int ndiel;  /**< Number of dielectric maps */
    char dielXpath[NOSH_MAXMOL][VMAX_ARGLEN];  /**< Paths to x-shifted
        dielectric map files */
    char dielYpath[NOSH_MAXMOL][VMAX_ARGLEN];  /**< Paths to y-shifted
        dielectric map files */
    char dielZpath[NOSH_MAXMOL][VMAX_ARGLEN];  /**< Paths to z-shifted
        dielectric map files */
    Vdata_Format dielfmt[NOSH_MAXMOL];  /**< Dielectric maps file formats */
    int nkappa;  /**< Number of kappa maps */
    char kappapath[NOSH_MAXMOL][VMAX_ARGLEN]; /**< Paths to kappa map files */
    Vdata_Format kappafmt[NOSH_MAXMOL];  /**< Kappa maps file formats */
    int npot;  /**< Number of potential maps */
    char potpath[NOSH_MAXMOL][VMAX_ARGLEN]; /**< Paths to potential map files */
    Vdata_Format potfmt[NOSH_MAXMOL];  /**< Potential maps file formats */
    int ncharge;  /**< Number of charge maps */
    char chargepath[NOSH_MAXMOL][VMAX_ARGLEN];   /**< Paths to charge map files */
    Vdata_Format chargefmt[NOSH_MAXMOL];  /**< Charge maps fileformats */
    int nmesh;  /**< Number of meshes */
    char meshpath[NOSH_MAXMOL][VMAX_ARGLEN];   /**< Paths to mesh files */
    Vdata_Format meshfmt[NOSH_MAXMOL];  /**< Mesh fileformats */
    int nprint;  /**< How many print sections? */
    NOsh_PrintType printwhat[NOSH_MAXPRINT];  /**< What do we print:  \li 0 =
        energy, \li 1 = force */
    int printnarg[NOSH_MAXPRINT];  /**< How many arguments in energy list */
    int printcalc[NOSH_MAXPRINT][NOSH_MAXPOP]; /**< ELEC id (see elec2calc) */
    int printop[NOSH_MAXPRINT][NOSH_MAXPOP];  /**< Operation id (0 = add, 1 =
        subtract) */
    int parsed;  /**< Have we parsed an input file yet? */
    char elecname[NOSH_MAXCALC][VMAX_ARGLEN]; /**< Optional user-specified name
        for ELEC statement */
    char apolname[NOSH_MAXCALC][VMAX_ARGLEN]; /**< Optional user-specified name
        for APOLAR statement */
};
```
##Multigrid Parameters
```C
struct sMGparm {

    MGparm_CalcType type;  /**< What type of MG calculation? */
    int parsed;  /**< Has this structure been filled? (0 = no, 1 = yes) */

    /* *** GENERIC PARAMETERS *** */
    int dime[3];  /**< Grid dimensions */
    int setdime;  /**< Flag, @see dime */
    Vchrg_Meth chgm;  /**< Charge discretization method */
    int setchgm;  /**< Flag, @see chgm */
    Vchrg_Src  chgs; /**< Charge source (Charge, Multipole, Induced Dipole,
                      * NL Induced */

    /* *** TYPE 0 PARAMETERS (SEQUENTIAL MANUAL) *** */
    int nlev;  /**< Levels in multigrid hierarchy
                *   @deprecated Just ignored now */
    int setnlev;  /**< Flag, @see nlev */
    double etol;  /**< User-defined error tolerance */
    int setetol;  /**< Flag, @see etol */
    double grid[3];  /**< Grid spacings */
    int setgrid;  /**< Flag, @see grid */
    double glen[3];  /**< Grid side lengths. */
    int setglen;  /**< Flag, @see glen */
    MGparm_CentMeth cmeth;  /**< Centering method */
    double center[3];  /**< Grid center. If ispart = 0, then this is
                        * only meaningful if cmeth = 0.  However, if
                        * ispart = 1 and cmeth = MCM_PNT, then this is the
                        * center of the non-disjoint (overlapping)
                        * partition.  If ispart = 1 and cmeth = MCM_MOL, then
                        * this is the vector that must be added to the
                        * center of the molecule to give the center of
                        * the non-disjoint partition.  */
    int centmol;  /**< Particular molecule on which we want to center the grid.
        This should be the appropriate index in an array of molecules, not the
        positive definite integer specified by the user. */
    int setgcent;  /**< Flag, @see cmeth */

    /* ******** TYPE 1 & 2 PARAMETERS (SEQUENTIAL & PARALLEL AUTO-FOCUS) *** */
    double cglen[3];  /**< Coarse grid side lengths */
    int setcglen;  /**< Flag, @see cglen */
    double fglen[3];  /**< Fine grid side lengths */
    int setfglen;  /**< Flag, @see fglen */
    MGparm_CentMeth ccmeth;  /**< Coarse grid centering method */
    double ccenter[3];  /**< Coarse grid center.  */
    int ccentmol;  /**< Particular molecule on which we want to center the grid.
        This should be the appropriate index in an array of molecules, not the
        positive definite integer specified by the user. */
    int setcgcent;  /**< Flag, @see ccmeth */
    MGparm_CentMeth fcmeth;  /**< Fine grid centering method */
    double fcenter[3];  /**< Fine grid center.  */
    int fcentmol; /**< Particular molecule on which we want to center the grid.
        This should be the appropriate index in an array of molecules, not the
        positive definite integer specified by the user. */
    int setfgcent;  /**< Flag, @see fcmeth */


    /* ********* TYPE 2 PARAMETERS (PARALLEL AUTO-FOCUS) ******** */
    double partDisjCenter[3];  /**< This gives the center
                                     of the disjoint partitions */
    double partDisjLength[3];  /**< This gives the lengths of the disjoint
                                * partitions */
    int partDisjOwnSide[6];  /**< Tells whether the boundary points are ours
                              * (1) or not (0) */

    int pdime[3];  /**< Grid of processors to be used in calculation */
    int setpdime;  /**< Flag, @see pdime */
    int proc_rank;  /**< Rank of this processor */
    int setrank;  /**< Flag, @see proc_rank */
    int proc_size;  /**< Total number of processors */
    int setsize;  /**< Flag, @see proc_size */
    double ofrac;  /**< Overlap fraction between procs */
    int setofrac;  /**< Flag, @see ofrac */
    int async; /**< Processor ID for asynchronous calculation */
    int setasync; /**< Flag, @see asynch */

    int nonlintype; /**< Linearity Type Method to be used */
    int setnonlintype; /**< Flag, @see nonlintype */

    int method;		/**< Solver Method */
    int setmethod; /**< Flag, @see method */

    int useAqua;  /**< Enable use of lpbe/aqua */
    int setUseAqua; /**< Flag, @see useAqua */
};
```
##PBE Parameters
```C
struct sPBEparm {

    int molid;  /**< Molecule ID to perform calculation on */
    int setmolid;  /**< Flag, @see molid */
    int useDielMap;  /**< Indicates whether we use external
                      * dielectric maps (note plural) */
    int dielMapID;  /**< Dielectric map ID (if used) */
    int useKappaMap;  /**< Indicates whether we use an external
                       * kappa map */
    int kappaMapID;  /**< Kappa map ID (if used) */
    int usePotMap;  /**< Indicates whether we use an external
                       * kappa map */
    int potMapID;  /**< Kappa map ID (if used) */

    int useChargeMap;  /**< Indicates whether we use an external
                        * charge distribution map */
    int chargeMapID;  /**< Charge distribution map ID (if used) */
    Vhal_PBEType pbetype;  /**< Which version of the PBE are we solving? */
    int setpbetype;  /**< Flag, @see pbetype */
    Vbcfl bcfl;  /**< Boundary condition method */
    int setbcfl;  /**< Flag, @see bcfl */
    int nion;  /**< Number of counterion species */
    int setnion;  /**< Flag, @see nion */
    double ionq[MAXION];  /**< Counterion charges (in e) */
    double ionc[MAXION];  /**< Counterion concentrations (in M) */
    double ionr[MAXION];  /**< Counterion radii (in A) */
    int setion[MAXION];  /**< Flag, @see ionq */
    double pdie;  /**< Solute dielectric */
    int setpdie;  /**< Flag, @see pdie */
    double sdens; /**< Vacc sphere density */
    int setsdens; /**< Flag, @see sdens */
    double sdie;  /**< Solvent dielectric */
    int setsdie;  /**< Flag, @see sdie */
    Vsurf_Meth srfm;  /**< Surface calculation method */
    int setsrfm;  /**< Flag, @see srfm */
    double srad;  /**< Solvent radius */
    int setsrad;  /**< Flag, @see srad */
    double swin;  /**< Cubic spline window */
    int setswin;  /**< Flag, @see swin */
    double temp;  /**< Temperature (in K) */
    int settemp;  /**< Flag, @see temp */

    double smsize; /**< SMPBE size */
    int setsmsize; /**< Flag, @see temp */

    double smvolume; /**< SMPBE size */
    int setsmvolume; /**< Flag, @see temp */

    PBEparm_calcEnergy calcenergy;  /**< Energy calculation flag */
    int setcalcenergy;  /**< Flag, @see calcenergy */
    PBEparm_calcForce calcforce;  /**< Atomic forces calculation */
    int setcalcforce;  /**< Flag, @see calcforce */

    /*----------------------------------------------------------------*/
    /* Added by Michael Grabe                                         */
    /*----------------------------------------------------------------*/

    double zmem;               /**< z value of membrane bottom */
    int setzmem;               /**< Flag */
    double Lmem;               /**< membrane width */
    int setLmem;               /**< Flag */
    double mdie;               /**< membrane dielectric constant */
    int setmdie;               /**< Flag */
    double memv;               /**< Membrane potential */
    int setmemv;               /**< Flag */

    /*----------------------------------------------------------------*/

    int numwrite;  /**< Number of write statements encountered */
    char writestem[PBEPARM_MAXWRITE][VMAX_ARGLEN]; /**< File stem to write
                                                    * data to */
    Vdata_Type writetype[PBEPARM_MAXWRITE];  /**< What data to write */
    Vdata_Format writefmt[PBEPARM_MAXWRITE];  /**< File format to write data
                                               * in */
    int writemat;  /**< Write out the operator matrix?
                    * \li 0 => no
                    * \li 1 => yes */
    int setwritemat;  /**< Flag, @see writemat */
    char writematstem[VMAX_ARGLEN];  /**< File stem to write mat */
    int writematflag;  /**< What matrix should we write:
                        * \li 0 => Poisson (differential operator)
                        * \li 1 => Poisson-Boltzmann operator linearized around
                        * solution (if applicable) */

    int parsed;  /**< Has this been filled with anything other
                  * than the default values? */
};
```
#Special Needs
* plugin to read dielectric maps
* plugin to read kappa maps
* plugin to read potential maps
* plugin to read charge maps

#Questions
* There is a check for a neutral charge in initmg that is commented out, but the calculation of the squared charge is still there.  D. Gohara.
* Aqua?
* Should output of remote jobs be routed back to invoker?  It seems to me a logical thing to do.
