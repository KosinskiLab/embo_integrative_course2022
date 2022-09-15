import IMP.pmi.io.crosslink
import IMP.pmi.macros
import IMP.pmi.output
import IMP.pmi.dof
import IMP.atom
import numpy


# Monkey patch CrossLinkDataBase.append_database to incorporate a post-2.11
# bug fix
def _append_database(self,CrossLinkDataBase2):
    name1=self.get_name()
    name2=CrossLinkDataBase2.get_name()
    if name1 == name2:
        name1=id(self)
        name2=id(CrossLinkDataBase2)
        self.set_name(name1)
        CrossLinkDataBase2.set_name(name2)
    #rename first database:
    new_data_base={}
    for k in self.data_base:
        new_data_base[k]=self.data_base[k]
    for k in CrossLinkDataBase2.data_base:
        new_data_base[k]=CrossLinkDataBase2.data_base[k]
    self.data_base=new_data_base
    self._update()


# Monkey patch AnalysisReplicaExchange.align to incorporate a post-2.11
# bug fix
def _align(self):
    tr = IMP.atom.get_transformation_aligning_first_to_second(self.sel1_alignment, self.sel0_alignment)
    for rb in self.rbs1:
        IMP.core.transform(rb, tr)
    for bead in self.beads1:
        if not IMP.core.NonRigidMember.get_is_setup(bead):
            IMP.core.transform(IMP.core.XYZ(bead), tr)  
    self.model.update()

# Monkey patch AnalysisReplicaExchange.__init__ to incorporate a post-2.11
# bug fix
def _init(self,model,
                 stat_files,
                 best_models=None,
                 score_key=None,
                 alignment=True):
    """
    Construction of the Class.
    @param model IMP.Model()
    @param stat_files list of string. Can be ascii stat files, rmf files names
    @param best_models Integer. Number of best scoring models, if None: all models will be read
    @param alignment boolean (Default=True). Align before computing the rmsd.
    """
    self.model=model
    self.best_models=best_models
    self.stath0=IMP.pmi.output.StatHierarchyHandler(model,stat_files,self.best_models,score_key,cache=True)
    self.stath1=IMP.pmi.output.StatHierarchyHandler(StatHierarchyHandler=self.stath0)

    self.rbs1, self.beads1 = IMP.pmi.tools.get_rbs_and_beads(IMP.pmi.tools.select_at_all_resolutions(self.stath1))
    self.rbs0, self.beads0 = IMP.pmi.tools.get_rbs_and_beads(IMP.pmi.tools.select_at_all_resolutions(self.stath0))
    self.sel0_rmsd=IMP.atom.Selection(self.stath0)
    self.sel1_rmsd=IMP.atom.Selection(self.stath1)
    self.sel0_alignment=IMP.atom.Selection(self.stath0)
    self.sel1_alignment=IMP.atom.Selection(self.stath1)
    self.clusters=[]
    # fill the cluster list with a single cluster containing all models
    c = IMP.pmi.output.Cluster(0)
    self.clusters.append(c)
    for n0,d0 in enumerate(self.stath0):
        c.add_member(n0,d0)
    self.pairwise_rmsd={}
    self.pairwise_molecular_assignment={}
    self.alignment=alignment
    self.symmetric_molecules={}
    self.issymmetricsel={}
    self.update_seldicts()
    self.molcopydict0=IMP.pmi.tools.get_molecules_dictionary_by_copy(IMP.atom.get_leaves(self.stath0))
    self.molcopydict1=IMP.pmi.tools.get_molecules_dictionary_by_copy(IMP.atom.get_leaves(self.stath1))


# This method isn't in IMP 2.11, so patch it in
def _classify_crosslinks_by_score(self, number_of_classes):
    '''Creates the requested number of classes and partitions crosslinks
       according to their identification scores. Classes are defined in
       the psi key.'''
    if self.id_score_key is not None:
        scores=self.get_values(self.id_score_key)
    else:
        raise ValueError('The crosslink database does not contain score values')
    minscore=min(scores)
    maxscore=max(scores)
    scoreclasses=numpy.linspace(minscore, maxscore, number_of_classes+1)
    if self.psi_key is None:
        self.create_new_keyword(self.psi_key,values_from_keyword=None)
    for xl in self:
        score=xl[self.id_score_key]
        for n,classmin in enumerate(scoreclasses[0:-1]):
            if score>=classmin and score<=scoreclasses[n+1]:
                xl[self.psi_key]=str("CLASS_"+str(n))
    self._update()

# This method isn't in IMP 2.11, so patch it in
def _get_output(self):
    return self[self.current_index].features

# This method isn't in IMP 2.11, so patch it in
def _write_seed(self, file_name, nreplicas):
    rmf_name = file_name
    o=IMP.pmi.output.Output()
    if nreplicas > len(self.stath1):
        raise ValueError('number of replicas exceeding the number of structures')
    o.init_rmf(rmf_name, [self.stath1], listofobjects=[self.stath1])

    structures=[]
    nloop=0
    while len(structures) < nreplicas:
        for c in self:
            try:
                structures.append(c.members[nloop])
            except:
                continue
        nloop+=1

    for n in structures:
        self.stath1[n]
        o.write_rmf(rmf_name)

    o.close_rmf(rmf_name)


# Monkey patch DegreesOfFreedom._setup_srb to incorporate a post-2.11 bug fix
def _setup_srb(self,hiers,max_trans,max_rot,axis):
    if axis is None:
        srbm = IMP.pmi.TransformMover(hiers[0][0].get_model(), max_trans, max_rot)
    else:
        srbm = IMP.pmi.TransformMover(hiers[0][0].get_model(),axis[0],axis[1],max_trans, max_rot)
    srbm.set_was_used(True)
    super_rigid_rbs,super_rigid_xyzs = IMP.pmi.tools.get_rbs_and_beads(hiers)
    ct = 0
    self.movers_particles_map[srbm]=[]
    for h in hiers:
        self.movers_particles_map[srbm]+=IMP.atom.get_leaves(h)
    for xyz in super_rigid_xyzs:
        srbm.add_xyz_particle(xyz)
        ct+=1
    for rb in super_rigid_rbs:
        srbm.add_rigid_body_particle(rb)
        ct+=1
    if ct>1:
        return srbm
    else:
        return 0

# Monkey patch DegreesOfFreedom._setup_srb to incorporate a post-2.11 bug fix
def _save_coordinates(self,cluster,rmf_name=None,reference="Absolute", prefix="./"):
    """
    Save the coordinates of the current cluster a single rmf file
    """
    print("saving coordinates",cluster)
    if self.alignment: self.set_reference(reference,cluster)
    o=IMP.pmi.output.Output()
    if rmf_name is None:
        rmf_name=prefix+'/'+str(cluster.cluster_id)+".rmf3"

    d1=self.stath1[cluster.members[0]]
    self.model.update()
    o.init_rmf(rmf_name, [self.stath1], listofobjects=[self.stath1])
    for n1 in cluster.members:
        d1=self.stath1[n1]
        self.model.update()
        self.apply_molecular_assignments(n1)
        if self.alignment: self.align()
        o.write_rmf(rmf_name)
        self.undo_apply_molecular_assignments(n1)
    o.close_rmf(rmf_name)


if IMP.__version__ == '2.11.1':
    IMP.pmi.macros.AnalysisReplicaExchange.__init__ = _init
    IMP.pmi.macros.AnalysisReplicaExchange.write_seed = _write_seed
    IMP.pmi.macros.AnalysisReplicaExchange.save_coordinates = _save_coordinates
    IMP.pmi.macros.AnalysisReplicaExchange.align = _align
    IMP.pmi.output.StatHierarchyHandler.get_output= _get_output
    IMP.pmi.io.crosslink.CrossLinkDataBase.append_database = _append_database
    IMP.pmi.io.crosslink.CrossLinkDataBase.classify_crosslinks_by_score \
        = _classify_crosslinks_by_score
    IMP.pmi.dof.DegreesOfFreedom._setup_srb = _setup_srb
