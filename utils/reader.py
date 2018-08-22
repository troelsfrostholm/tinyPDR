import numpy as np

class Run:
  """ Encapsulates output of a tinyPDR run """

  def __init__(self, filename):
    with open(filename, 'r') as f:
      header = f.readline().split()
      self.ngrid, self.ntime, self.nmols, self.nbins = map(int, header)
      self.data = np.loadtxt(f)
      self.offset_mol = 5
      self.offset_tau = self.offset_mol + self.nmols
      self.data.shape=(self.ntime, self.ngrid, self.offset_mol+self.nmols+self.nbins)

  def getSnapshot(self, itime):
    return self.data[itime,:,:]

  def getMolecules(self, itime=None):
    mols = self.data[:,:,self.offset_mol:self.offset_mol+self.nmols]
    if itime:
      return mols[itime,:,:]
    return mols

  def getTau(self, itime=None):
    tau = self.data[:,:,self.offset_tau:]
    if itime:
      return tau[itime,:,:]
    return tau

  def getTime(self, itime=None):
    time = self.data[:,:,0]
    if itime:
      return time[itime,:]
    return time

  def getR(self, itime=None):
    r = self.data[:,:,1]
    if itime:
      return r[itime,:]
    return r

  def getAv(self, itime=None):
    Av = self.data[:,:,2]
    if itime:
      return Av[itime,:]
    return Av

  def getTgas(self, itime=None):
    Tgas = self.data[:,:,3]
    if itime:
      return Tgas[itime,:]
    return Tgas

  def getTdust(self, itime=None):
    Tdust = self.data[:,:,4]
    if itime:
      return Tdust[itime,:]
    return Tdust

class FluxFile:
  def __init__(self, filename):
    with open(filename, 'r') as f:
      header = f.readline().split()
      self.ngrid, self.ntime, self.nrea = map(int, header)
      self.data = np.loadtxt(f)
      self.offset_rea = 3
      self.data.shape=(self.ntime, self.ngrid, self.offset_rea+self.nrea)

    self.read_reaction_verbatim()
    self.time = self.data[:,0,0]
    self.distance = self.data[0,:,1]
    self.av = self.data[0,:,2]
    self.reactions = Reactions(self.reaction_names, self.data[:,:,self.offset_rea:])

  def getReactions(self):
    return self.reactions

  def read_reaction_verbatim(self, filename="reactions_verbatim.dat"):
    with open(filename, 'r') as f:
      self.reaction_names = np.array(f.read().splitlines())

class Reactions:

  def __init__(self,reaction_names, fluxes):
    self.reaction_names = reaction_names
    self.fluxes = fluxes

  def selectByMolecule(self, molecule):
    ids = [id for (id, name) 
              in zip(range(len(self.reaction_names)), self.reaction_names)
              if Reactions.reactionContainsMolecule(name, molecule)]
    return Reactions(self.reaction_names[ids], self.fluxes[:,:,ids])

  def selectByTimeStep(self, ll, ul):
    return Reactions(self.reaction_names, self.fluxes[ll:ul,:,:])

  def selectByCell(self, ll, ul):
    return Reactions(self.reaction_names, self.fluxes[:,ll:ul,:])

  def top(self, N):
    itop = frozenset()
    s = self.fluxes.shape
    for itime in xrange(s[0]):
      for igrid in xrange(s[1]):
        fluxes = self.fluxes[itime,igrid,:]
        i = np.argsort(-fluxes)
        itop = itop | frozenset(i[:N])
    itop = list(itop)
    return Reactions(self.reaction_names[itop], self.fluxes[:,:,itop])

  def sorted(self):
    i = np.argsort(-self.fluxes)
    return Reactions(self.reaction_names[i], self.fluxes[:,:,i])

  def iter(self):
    current = 0
    while current < len(self.reaction_names):
      yield (self.reaction_names[current], self.fluxes[:,:,current])
      current += 1

  @staticmethod
  def reactionContainsMolecule(reaction_name, molecule):
    import re
    molecule = re.escape(molecule)
    p1 = re.search('^'+molecule+' ', reaction_name)
    p2 = re.search(' '+molecule+'$', reaction_name)
    p3 = re.search(' '+molecule+' ', reaction_name)
    if p1 or p2 or p3:
      return True
    return False