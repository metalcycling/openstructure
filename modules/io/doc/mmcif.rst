mmCIF File Format
--------------------------------------------------------------------------------

.. currentmodule:: ost.io

The mmCIF file format is a container for structural entities provided by the
PDB. Here we describe how to load those files and how to deal with information
provided above the legacy PDB format (:class:`MMCifInfo`,
:class:`MMCifInfoCitation`, :class:`MMCifInfoTransOp`,
:class:`MMCifInfoBioUnit`, :class:`MMCifInfoStructDetails`,
:class:`MMCifInfoObsolete`, :class:`MMCifInfoStructRef`,
:class:`MMCifInfoStructRefSeq`, :class:`MMCifInfoStructRefSeqDif`,
:class:`MMCifInfoRevisions`, :class:`MMCifInfoEntityBranchLink`).


Loading mmCIF Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: ost.io.LoadMMCIF


Categories Available
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following categories of a mmCIF file are considered by the reader:

* ``atom_site``: Used to build the :class:`~ost.mol.EntityHandle`
* ``entity``: Involved in setting :class:`~ost.mol.ChainType` of chains
* ``entity_poly``: Involved in setting :class:`~ost.mol.ChainType` of chains
* ``citation``: Goes into :class:`MMCifInfoCitation`
* ``citation_author``: Goes into :class:`MMCifInfoCitation`
* ``exptl``: Goes into :class:`MMCifInfo` as :attr:`~MMCifInfo.method`.
* ``refine``: Goes into :class:`MMCifInfo` as :attr:`~MMCifInfo.resolution`,
  :attr:`~MMCifInfo.r_free` and :attr:`~MMCifInfo.r_work`.
* ``pdbx_struct_assembly``: Used for :class:`MMCifInfoBioUnit`.
* ``pdbx_struct_assembly_gen``: Used for :class:`MMCifInfoBioUnit`.
* ``pdbx_struct_oper_list``: Used for :class:`MMCifInfoBioUnit`.
* ``struct``: Details about a structure, stored in
  :class:`MMCifInfoStructDetails`.
* ``struct_conf``: Stores secondary structure information (practically helices)
  in the :class:`~ost.mol.EntityHandle`
* ``struct_sheet_range``: Stores secondary structure information for sheets in
  the :class:`~ost.mol.EntityHandle`
* ``pdbx_database_PDB_obs_spr``: Verbose information on obsoleted/ superseded
  entries, stored in :class:`MMCifInfoObsolete`
* ``struct_ref`` stored in :class:`MMCifInfoStructRef`
* ``struct_ref_seq`` stored in :class:`MMCifInfoStructRefSeqDif`
* ``struct_ref_seq_dif`` stored in :class:`MMCifInfoStructRefDif`
* ``database_pdb_rev`` (mmCIF dictionary version < 5) stored in
  :class:`MMCifInfoRevisions`
* ``pdbx_audit_revision_history`` and ``pdbx_audit_revision_details``
  (mmCIF dictionary version >= 5) used to fill :class:`MMCifInfoRevisions`
* ``pdbx_entity_branch`` and ``pdbx_entity_branch_link`` used for
  :class:`MMCifInfoEntityBranchLink`, a list of links is available by
  :meth:`~MMCifInfo.GetEntityBranchLinks` and
  :meth:`~MMCifInfo.GetEntityBranchByChain`

Notes:

* Structures in mmCIF format can have two chain names. The "new" chain name
  extracted from ``atom_site.label_asym_id`` is used to name the chains in the
  :class:`~ost.mol.EntityHandle`. The "old" (author provided) chain name is
  extracted from ``atom_site.auth_asym_id`` for the first atom of the chain.
  It is added as string property named "pdb_auth_chain_name" to the
  :class:`~ost.mol.ChainHandle`. The mapping is also stored in
  :class:`MMCifInfo` as :meth:`~MMCifInfo.GetMMCifPDBChainTr` and
  :meth:`~MMCifInfo.GetPDBMMCifChainTr` if SEQRES records are read in
  :func:`~ost.io.LoadMMCIF` and a non-empty SEQRES record exists for that chain
  (this should exclude ligands and water).
* Molecular entities in mmCIF are identified by an ``entity.id``. Each chain is
  mapped to an ID in :class:`MMCifInfo` as :meth:`~MMCifInfo.GetMMCifEntityIdTr`.

Info Classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Information from mmCIF files that goes beyond structural data, is kept in a
special container, the :class:`MMCifInfo` class. Here is a detailed description
of the annotation available.

.. class:: MMCifInfo

  This is the container for all bits of non-molecular data pulled from a mmCIF
  file.

  .. attribute:: citations

    Stores a list of citations (:class:`MMCifInfoCitation`).

    Also available as :meth:`GetCitations`.

  .. attribute:: biounits

    Stores a list of biounits (:class:`MMCifInfoBioUnit`).

    Also available as :meth:`GetBioUnits`.

  .. attribute:: method

    Stores the experimental method used to create the structure.

    Also available as :meth:`GetMethod`. May also be modified by
    :meth:`SetMethod`.

  .. attribute:: resolution

    Stores the resolution of the crystal structure. Set to 0 if no value in
    loaded mmCIF file.

    Also available as :meth:`GetResolution`. May also be modified by
    :meth:`SetResolution`.

  .. attribute:: r_free

    Stores the R-free value of the crystal structure. Set to 0 if no value in
    loaded mmCIF file.

    Also available as :meth:`GetRFree`. May also be modified by
    :meth:`SetRFree`.

  .. attribute:: r_work

    Stores the R-work value of the crystal structure. Set to 0 if no value in
    loaded mmCIF file.

    Also available as :meth:`GetRWork`. May also be modified by
    :meth:`SetRWork`.

  .. attribute:: operations

    Stores the operations needed to transform a crystal structure into a
    bio unit.

    Also available as :meth:`GetOperations`. May also be modified by
    :meth:`AddOperation`.

  .. attribute:: struct_details

    Stores details about the structure in a :class:`MMCifInfoStructDetails`
    object.

    Also available as :meth:`GetStructDetails`. May also be modified by
    :meth:`SetStructDetails`.

  .. attribute:: struct_refs

    Lists all links to external databases in the mmCIF file.

  .. attribute:: revisions

    Stores a simple history of a PDB entry.

    Also available as :meth:`GetRevisions`. May be extended by
    :meth:`AddRevision`.

    :type: :class:`MMCifInfoRevisions`

  .. attribute:: obsolete

    Stores information about obsoleted / superseded entries.

    Also available as :meth:`GetObsoleteInfo`. May also be modified by
    :meth:`SetObsoleteInfo`.

    :type: :class:`MMCifInfoObsolete`

  .. method:: AddCitation(citation)

    Add a citation to the citation list of an info object.

    :param citation: Citation to be added.
    :type citation: :class:`MMCifInfoCitation`


  .. method:: AddAuthorsToCitation(id, authors)

    Adds a list of authors to a specific citation.

    :param id: Identifier of the citation.
    :type id: :class:`str`
    :param authors: List of authors.
    :type authors: :class:`~ost.StringList`

  .. method:: GetCitations()

    See :attr:`citations`

  .. method:: AddBioUnit(biounit)

    Add a bio unit to the bio unit list of an info object. If the
    :attr:`id <MMCifInfoBioUnit.id>` of ``biounit`` already exists in the set
    of assemblies, both will be merged. This means that
    :attr:`chain <MMCifInfoBioUnit.chains>` and
    :attr:`operations <MMCifInfoBioUnit.operations>` lists will be concatenated
    and the interval lists
    (:attr:`operationsintervalls <MMCifInfoBioUnit.operationsintervalls>`,
    :attr:`chainintervalls <MMCifInfoBioUnit.chainintervalls>`) will be updated.

    :param biounit: Bio unit to be added.
    :type biounit: :class:`MMCifInfoBioUnit`

  .. method:: GetBioUnits()

    See :attr:`biounits`

  .. method:: SetMethod(method)

    See :attr:`method`

  .. method:: GetMethod()

    See :attr:`method`

  .. method:: SetResolution(resolution)

    See :attr:`resolution`

  .. method:: GetResolution()

    See :attr:`resolution`

  .. method:: AddOperation(operation)

    See :attr:`operations`

  .. method:: GetOperations()

    See :attr:`operations`

  .. method:: SetStructDetails(details)

    See :attr:`struct_details`

  .. method:: GetStructDetails()

  .. method:: AddMMCifPDBChainTr(cif_chain_id, pdb_chain_id)

    Set up a translation for a certain mmCIF chain name to the traditional PDB
    chain name.

    :param cif_chain_id: atom_site.label_asym_id
    :type cif_chain_id: :class:`str`
    :param pdb_chain_id: atom_site.auth_asym_id
    :type pdb_chain_id: :class:`str`

  .. method:: GetMMCifPDBChainTr(cif_chain_id)

    Get the translation of a certain mmCIF chain name to the traditional PDB
    chain name.

    :param cif_chain_id: atom_site.label_asym_id
    :type cif_chain_id: :class:`str`
    :returns: atom_site.auth_asym_id as :class:`str` (empty if no mapping)

  .. method:: AddPDBMMCifChainTr(pdb_chain_id, cif_chain_id)

    Set up a translation for a certain PDB chain name to the mmCIF chain name.

    :param pdb_chain_id: atom_site.auth_asym_id
    :type pdb_chain_id: :class:`str`
    :param cif_chain_id: atom_site.label_asym_id
    :type cif_chain_id: :class:`str`

  .. method:: GetPDBMMCifChainTr(pdb_chain_id)

    Get the translation of a certain PDB chain name to the mmCIF chain name.

    :param pdb_chain_id: atom_site.auth_asym_id
    :type pdb_chain_id: :class:`str`
    :returns: atom_site.label_asym_id as :class:`str` (empty if no mapping)

  .. method:: AddMMCifEntityIdTr(cif_chain_id, entity_id)

    Set up a translation for a certain mmCIF chain name to the mmCIF entity ID.

    :param cif_chain_id: atom_site.label_asym_id
    :type cif_chain_id: :class:`str`
    :param entity_id: atom_site.label_entity_id
    :type entity_id: :class:`str`

  .. method:: GetMMCifEntityIdTr(cif_chain_id)

    Get the translation of a certain mmCIF chain name to the mmCIF entity ID.

    :param cif_chain_id: atom_site.label_asym_id
    :type cif_chain_id: :class:`str`
    :returns: atom_site.label_entity_id as :class:`str` (empty if no mapping)

  .. method:: AddRevision(num, date, status, major=-1, minor=-1)

    Add a new iteration to the revision history.
    See :meth:`MMCifInfoRevisions.AddRevision`.

  .. method:: GetRevisions()

    See :attr:`revisions`

  .. method:: SetRevisionsDateOriginal(date)

    Set the date, when this entry first entered the PDB. Ignored if it was set
    in the past. See :meth:`MMCifInfoRevisions.SetDateOriginal`.

  .. method:: GetObsoleteInfo()

    See :attr:`obsolete`

  .. method:: SetObsoleteInfo()

    See :attr:`obsolete`

  .. method:: GetEntityBranchLinks()

    Get bond information for branched entities. Returns all
    :class:`MMCifInfoEntityBranchLink` objects in one list. Chain and residue
    information is available by the stored
    :class:`AtomHandles <ost.mol.AtomHandle>` of each entry.

    :returns: :class:`list` of :class:`MMCifInfoEntityBranchLink`

  .. method:: GetEntityBranchByChain(chain_name)

    Get bond information for chains with branched entities. Returns all
    :class:`MMCifInfoEntityBranchLink` objects in one list if chain is a
    branched entity, an empty list otherwise.

    :param chain_name: Chain name to check for branch links
    :type chain_name: :class:`str`
    :returns: :class:`list` of :class:`MMCifInfoEntityBranchLink`

  .. method:: AddEntityBranchLink(chain_name, atom1, atom2, bond_order)

    Add bond information for a branched entity.

    :param chain_name: Chain the bond belongs to
    :type chain_name: :class:`str`
    :param atom1: First atom of the bond
    :type atom1: :class:`~ost.mol.AtomHandle`
    :param atom2: Second atom of the bond
    :type atom2: :class:`~ost.mol.AtomHandle`
    :param bond_order: Bond order (e.g. 1=single, 2=double, 3=triple)
    :type bond_order: :class:`int`
    :returns: Nothing

  .. method:: GetEntityBranchChainNames

    Get a list of chain names which contain branched entities.

    :returns: :class:`list` of :class:`str`

  .. method:: GetEntityBranchChains

    Get a list of chains which contain branched entities.

    :returns: :class:`list` of :class:`~ost.mol.ChainHandle`

  .. method:: ConnectBranchLinks

    Establish all bonds stored for branched entities.

.. class:: MMCifInfoCitation

  This stores citation information from an input file.

  .. attribute:: id

    Stores an internal identifier for a citation. If not provided, resembles an
    empty string.

    Also available as :meth:`GetID`. May also be modified by :meth:`SetID`.

  .. attribute:: cas

    Stores a Chemical Abstract Service identifier if available. If not
    provided, resembles an empty string.

    Also available as :meth:`GetCAS`. May also be modified by :meth:`SetCas`.

  .. attribute:: isbn

    Stores the ISBN code, presumably for cited books.  If not
    provided, resembles an empty string.

    Also available as :meth:`GetISBN`. May also be modified by :meth:`SetISBN`.

  .. attribute:: published_in

    Stores the book or journal title of a publication. Should take the full
    title, no abbreviations. If not provided, resembles an empty string.

    Also available as :meth:`GetPublishedIn`. May also be modified by
    :meth:`SetPublishedIn`.

  .. attribute:: volume

    Supposed to store volume information for journals. Since the volume number
    is not always a simple integer, it is stored as a string. If not provided,
    resembles an empty string.

    Also available as :meth:`GetVolume`. May also be modified by
    :meth:`SetVolume`.

  .. attribute:: page_first

    Stores the first page of a publication. Since the page numbers are not
    always a simple integers, they are stored as strings. If not provided,
    resembles empty strings.

    Also available as :meth:`GetPageFirst`. May also be modified by
    :meth:`SetPageFirst`.

  .. attribute:: page_last

    Stores the last page of a publication. Since the page numbers are not
    always a simple integers, they are stored as strings. If not provided,
    resembles empty strings.

    Also available as :meth:`GetPageLast`. May also be modified by
    :meth:`SetPageLast`.

  .. attribute:: doi

    Stores the Document Object Identifier as used by doi.org for a cited
    document. If not provided, resembles an empty string.

    Also available as :meth:`GetDOI`. May also be modified by :meth:`SetDOI`.

  .. attribute:: pubmed

    Stores the PubMed accession number. If not provided, is set to 0.

    Also available as :meth:`GetPubMed`. May also be modified by
    :meth:`SetPubmed`.

  .. attribute:: year

    Stores the publication year. If not provided, is set to 0.

    Also available as :meth:`GetYear`. May also be modified by :meth:`SetYear`.

  .. attribute:: title

    Stores a title. If not provided, is set to an empty string.

    Also available as :meth:`GetTitle`. May also be modified by
    :meth:`SetTitle`.

  .. attribute:: book_publisher

    Name of publisher of the citation, relevant for books and book chapters.

    Also available as :meth:`GetBookPublisher` and :meth:`SetBookPublisher`.

  .. attribute:: book_publisher_city

    City of the publisher of the citation, relevant for books and book
    chapters.

    Also available as :meth:`GetBookPublisherCity` and
    :meth:`SetBookPublisherCity`.

 
  .. attribute:: citation_type

     Defines where a citation was published. Either journal, book or unknown.

     Also available as :meth:`GetCitationType`. May also be modified by
     :meth:`SetCitationType` with values from :class:`MMCifInfoCType`. For
     conveinience setters :meth:`SetCitationTypeJournal`,
     :meth:`SetCitationTypeBook` and :meth:`SetCitationTypeUnknown` exist.

     For checking the type of a citation, :meth:`IsCitationTypeJournal`,
     :meth:`IsCitationTypeBook` and :meth:`IsCitationTypeUnknown` can be used.

  .. attribute:: authors

    Stores a :class:`~ost.StringList` of authors.

    Also available as :meth:`GetAuthorList`. May also be modified by
    :meth:`SetAuthorList`.

  .. method:: GetCAS()
    
    See :attr:`cas`

  .. method:: SetCAS(cas)

    See :attr:`cas`

  .. method:: GetISBN()
    
    See :attr:`isbn`

  .. method:: SetISBN(isbn)

    See :attr:`isbn`

  .. method:: GetPublishedIn()
    
    See :attr:`published_in`

  .. method:: SetPublishedIn(title)

    See :attr:`published_in`

  .. method:: GetVolume()
    
    See :attr:`volume`

  .. method:: SetVolume(volume)

    See :attr:`volume`

  .. method:: GetPageFirst()
    
    See :attr:`page_first`

  .. method:: SetPageFirst(first)

    See :attr:`page_first`

  .. method:: GetPageLast()
    
    See :attr:`page_last`

  .. method:: SetPageLast(last)

    See :attr:`page_last`

  .. method:: GetDOI()
    
    See :attr:`doi`

  .. method:: SetDOI(doi)

    See :attr:`doi`

  .. method:: GetPubMed()
    
    See :attr:`pubmed`

  .. method:: SetPubMed(no)

    See :attr:`pubmed`

  .. method:: GetYear()
    
    See :attr:`year`

  .. method:: SetYear(year)

    See :attr:`year`

  .. method:: GetTitle()
    
    See :attr:`title`

  .. method:: SetTitle(title)

    See :attr:`title`

  .. method:: GetBookPublisher

    See :attr:`book_publisher`

  .. method:: SetBookPublisher

    See :attr:`book_publisher`

  .. method:: GetBookPublisherCity

    See :attr:`book_publisher_city`

  .. method:: SetBookPublisherCity

    See :attr:`book_publisher_city`

  .. method:: GetCitationType()

    See :attr:`citation_type`

  .. method:: SetCitationType(publication_type)

    See :attr:`citation_type`

  .. method:: SetCitationTypeJournal()

    See :attr:`citation_type`

  .. method:: SetCitationTypeBook()

    See :attr:`citation_type`

  .. method:: SetCitationTypeUnknown()

    See :attr:`citation_type`

  .. method:: IsCitationTypeJournal()

    See :attr:`citation_type`

  .. method:: IsCitationTypeBook()

    See :attr:`citation_type`

  .. method:: IsCitationTypeUnknown()

    See :attr:`citation_type`
 
  .. method:: GetAuthorList()

    See :attr:`authors`

  .. method:: SetAuthorList(list)

    See :attr:`authors`


.. class:: MMCifInfoTransOp

  This stores operations needed to transform an :class:`~ost.mol.EntityHandle`
  into a bio unit.

  .. attribute:: id

    A unique identifier. If not provided, resembles an empty string.

    Also available as :meth:`GetID`. May also be modified by
    :meth:`SetID`.

  .. attribute:: type

    Describes the operation. If not provided, resembles an empty string.

    Also available as :meth:`GetType`. May also be modified by
    :meth:`SetType`.

  .. attribute:: translation

    The translational vector. Also available as :meth:`GetVector`. May also be

    modified by :meth:`SetVector`.

  .. attribute:: rotation

    The rotational matrix. Also available as :meth:`GetMatrix`. May also be

    modified by :meth:`SetMatrix`.

  .. method:: GetID()

    See :attr:`id`

  .. method:: SetID(id)

    See :attr:`id`

  .. method:: GetType()

    See :attr:`type`

  .. method:: SetType(type)

    See :attr:`type`

  .. method:: GetVector()

    See :attr:`translation`

  .. method:: SetVector(x, y, z)

    See :attr:`translation`

  .. method:: GetMatrix()

    See :attr:`rotation`

  .. method:: SetMatrix(i00,i01, i02, i10,i11, i12, i20,i21, i22)

    See :attr:`rotation`

.. class:: MMCifInfoBioUnit

  This stores information how a structure is to be assembled to form the
  bio unit.

  .. attribute:: id

    The id of a bio unit as given by the original mmCIF file.

    Also available as :meth:`GetID`. May also be modified by :meth:`SetID`.

    :type: :class:`str`

  .. attribute:: details

    Special aspects of the biological assembly. If not provided, resembles an
    empty string.

    Also available as :meth:`GetDetails`. May also be modified by
    :meth:`SetDetails`.

  .. attribute:: method_details

    Details about the method used to determine this biological assembly.

    Also available as :meth:`GetMethodDetails`. May also be modified by
    :meth:`SetMethodDetails`.

  .. attribute:: chains

    Chains involved in this bio unit. If not provided, resembles an empty list.

    Also available as :meth:`GetChainList`. May also be modified by
    :meth:`AddChain` or :meth:`SetChainList`.

  .. attribute:: chainintervals

    List of intervals on the chain list. Needed if there a several sets of
    chains and transformations to create the bio unit. Comes as a list of
    tuples. First component is the start, second is the right border of the
    interval.

    Also available as :meth:`GetChainIntervalList`. Is automatically modified by
    :meth:`AddChain`, :meth:`SetChainList` and :meth:`MMCifInfo.AddBioUnit`.

  .. attribute:: operations

    Translations and rotations needed to create the bio unit. Filled with
    objects of class :class:`MMCifInfoTransOp`.

    Also available as :meth:`GetOperations`. May be modified by
    :meth:`AddOperations`

  .. attribute:: operationsintervalls

    List of intervals on the operations list. Needed if there a several sets of
    chains and transformations to create the bio unit. Comes as a list of
    tuples. First component is the start, second is the right border of the
    interval.

    Also available as :meth:`GetOperationsIntervalList`. Is automatically
    modified by :meth:`AddOperations` and :meth:`MMCifInfo.AddBioUnit`.

  .. method:: GetID()

    See :attr:`id`

  .. method:: SetID(id)

    See :attr:`id`

  .. method:: GetDetails()

    See :attr:`details`

  .. method:: SetDetails(details)

    See :attr:`details`

  .. method:: GetMethodDetails()

    See :attr:`method_details`

  .. method:: SetMethodDetails(details)

    See :attr:`method_details`

  .. method:: GetChainList()

    See :attr:`chains`

  .. method:: SetChainList(chains)

    See :attr:`chains`, also resets :attr:`chainintervalls` to contain only one
    interval enclosing the whole chain list.

    :param chains: List of chain names.
    :type chains: :class:`~ost.StringList`

  .. method:: AddChain(chain name)

    See :attr:`chains`, also extends the right border of the last entry in
    :attr:`chainintervalls`.

  .. method:: GetChainIntervalList()

    See :attr:`chainintervals`

  .. method:: GetOperations()

    See :attr:`operations`

  .. method:: AddOperations(list of operations)

    See :attr:`operations`, also extends the right border of the last entry in
    :attr:`operationsintervalls`.

  .. method:: GetOperationsIntervalList()

    See :attr:`operationsintervalls`

  .. function:: PDBize(asu, seqres=None, min_polymer_size=None, transformation=False, peptide_min_size=10, nucleicacid_min_size=10, saccharide_min_size=10)

    Returns the biological assembly (bio unit) for an entity. The new entity
    created is well suited to be saved as a PDB file. Therefore the function
    tries to meet the requirements of single-character chain names. The
    following measures are taken.
  
      - All ligands are put into one chain (_)
      - Water is put into one chain (-)
      - Each polymer gets its own chain, named A-Z 0-9 a-z.
      - The description of non-polymer chains will be put into a generic string
        property called description on the residue level.
      - Ligands that resemble a polymer but have less than *min_polymer_size* /
        *peptide_min_size* / *nucleicacid_min_size* / *saccharide_min_size*
        residues are assigned the same numeric residue number. The residues are
        distinguished by insertion code.
      - Sometimes bio units exceed the coordinate system storable in a PDB file.
        In that case, the box around the entity will be aligned to the lower
        left corner of the coordinate system.

    Since this function is at the moment mainly used to create biounits from
    mmCIF files to be saved as PDBs, the function assumes that the
    :class:`~ost.mol.ChainType` properties are set correctly. 

    :param asu:  Asymmetric unit to work on. Should be created from a mmCIF
                 file.
    :type asu: :class:`~ost.mol.EntityHandle`
    :param seqres: If set to a valid sequence list, the length of the seqres 
      records will be used to determine if a certain chain has the minimally 
      required length.
    :type seqres: :class:`~ost.seq.SequenceList`
    :param min_polymer_size:  The minimal number of residues a polymer needs to 
      get its own chain. Everything below that number will be sorted into the 
      ligand chain. Overrides *peptide_min_size*, *nucleicacid_min_size* and
      *saccharide_min_size* if set to a value different than None.
    :type min_polymer_size: int
    :param transformation:  If set, return the transformation matrix used to
      move the bounding box of the bio unit to the lower left corner.
    :type transformation: :class:`bool`
    :param peptide_min_size: Minimal size to get an individual chain for a
      polypeptide. Is overridden by *min_polymer_size*.
    :type peptide_min_size: :class:`int`
    :param nucleicacid_min_size: Minimal size to get an individual chain for a
      polynucleotide. Is overridden by *min_polymer_size*.
    :type nucleicacid_min_size: :class:`int`
    :param saccharide_min_size: Minimal size to get an individual chain for an
      oligosaccharide or polysaccharide. Is overridden by *min_polymer_size*.
    :type saccharide_min_size: :class:`int`

.. class:: MMCifInfoStructDetails

  Holds details about the structure.

  .. attribute:: entry_id

    Identifier for a curtain data block. If not provided, resembles an empty
    string.

    Also available as :meth:`GetEntryID`. May also be modified by
    :meth:`SetEntryID`.

  .. attribute:: title

    Set a title for the structure.

    Also available as :meth:`GetTitle`. May also be modified by
    :meth:`SetTitle`.

  .. attribute:: casp_flag

    Tells whether this structure was a target in some competition.

    Also available as :meth:`GetCASPFlag`. May also be modified by
    :meth:`SetCASPFlag`.

  .. attribute:: descriptor

    Descriptor for an NDB structure or the unstructured content of a PDB COMPND
    record.

    Also available as :meth:`GetDescriptor`. May also be modified by
    :meth:`SetDescriptor`.

  .. attribute:: mass

    Molecular mass of a molecule.

    Also available as :meth:`GetMass`. May also be modified by
    :meth:`SetMass`.

  .. attribute:: mass_method

    Method used to determine the molecular weight.

    Also available as :meth:`GetMassMethod`. May also be modified by
    :meth:`SetMassMethod`.

  .. attribute:: model_details

    Details about how the structure was determined.

    Also available as :meth:`GetModelDetails`. May also be modified by
    :meth:`SetModelDetails`.

  .. attribute:: model_type_details

    Details about how the type of the structure was determined.

    Also available as :meth:`GetModelTypeDetails`. May also be modified by
    :meth:`SetModelTypeDetails`.

  .. method:: GetEntryID()

    See :attr:`entry_id`

  .. method:: SetEntryID(id)

    See :attr:`entry_id`

  .. method:: GetTitle()

    See :attr:`title`

  .. method:: SetTitle(title)

    See :attr:`title`

  .. method:: GetCASPFlag()

    See :attr:`casp_flag`

  .. method:: SetCASPFlag(flag)

    See :attr:`casp_flag`

  .. method:: GetDescriptor()

    See :attr:`descriptor`

  .. method:: SetDescriptor(descriptor)

    See :attr:`descriptor`

  .. method:: GetMass()

    See :attr:`mass`

  .. method:: SetMass(mass)

    See :attr:`mass`

  .. method:: GetMassMethod()

    See :attr:`mass_method`

  .. method:: SetMassMethod(method)

    See :attr:`mass_method`

  .. method:: GetModelDetails()

    See :attr:`model_details`

  .. method:: SetModelDetails(details)

    See :attr:`model_details`

  .. method:: GetModelTypeDetails()

    See :attr:`model_type_details`

  .. method:: SetModelTypeDetails(details)

    See :attr:`model_type_details`

.. class:: MMCifInfoObsolete

  Holds details on obsolete / superseded structures. The data is
   available both in the obsolete and in the replacement entries.

  .. attribute:: date

    When was the entry replaced?

    Also available as :meth:`GetDate`. May also be modified by
    :meth:`SetDate`.

  .. attribute:: id

    Type of change. Either *Obsolete* or *Supersede*. Returns a string starting
    upper case. Has to be set via ``OBSLTE`` or ``SPRSDE``.

    Also available as :meth:`GetID`. May also be modified by
    :meth:`SetID`.

  .. attribute:: pdb_id

    ID of the replacing entry.

    Also available as :meth:`GetPDBID`. May also be modified by
    :meth:`SetPDBID`.

  .. attribute:: replace_pdb_id

    ID of the replaced entry.

    Also available as :meth:`GetReplacedPDBID`. May also be modified by
    :meth:`SetReplacedPDBID`.

  .. method:: GetDate()

    See :attr:`date`

  .. method:: SetDate(date)

    See :attr:`date`

  .. method:: GetID()

    See :attr:`id`

  .. method:: SetID(id)

    See :attr:`id`

  .. method:: GetPDBID()

    See :attr:`pdb_id`

  .. method:: SetPDBID(flag)

    See :attr:`pdb_id`

  .. method:: GetReplacedPDBID()

    See :attr:`replace_pdb_id`

  .. method:: SetReplacedPDBID(descriptor)

    See :attr:`replace_pdb_id`

.. class:: MMCifInfoStructRef

  Holds the information of the struct_ref category. The category describes the 
  link of polymers in the mmCIF file to sequences stored in external databases 
  such as UniProt. The related categories ``struct_ref_seq`` and 
  ``struct_ref_seq_dif`` also list differences between the sequences of the 
  deposited structure and the sequences in the database. Two prominent examples 
  of such differences are point mutations and/or expression tags.

  .. attribute:: db_name

    
    Name of the external database, for example UNP for UniProt.

    :type: :class:`str`


  .. attribute:: db_id
    
    Name of the reference sequence in the database pointed to by :attr:`db_name`.

    :type: :class:`str`
  
  .. attribute:: db_access
    
    Alternative accession code for the sequence in the database pointed to by 
    :attr:`db_name`.

    :type: :class:`str`

  .. method:: GetAlignedSeq(name)

    Returns the aligned sequence for the given name, None if the sequence does 
    not exist.
  
  .. attribute:: aligned_seqs

    List of aligned sequences (all entries of the struct_ref_seq category 
    mapping to this struct_ref).

.. class:: MMCifInfoStructRefSeq

  An aligned range of residues between a sequence in a reference database and the 
  deposited sequence.

  .. attribute:: align_id
    
    Uniquely identifies every struct_ref_seq item in the mmCIF file.

    :type: :class:`str`

  .. attribute:: seq_begin
                 seq_end

    The starting point (1-based) and end point of the aligned range in the 
    deposited sequence, respectively.

    :type: :class:`int`
   
  .. attribute:: db_begin
                 db_end

    The starting point (1-based) and end point of the aligned range in the 
    database sequence, respectively.

    :type: :class:`int`

  .. attribute:: difs

    List of differences between the deposited sequence and the sequence in the 
    database.

  .. attribute:: chain_name

    Chain name of the polymer in the mmCIF file.

.. class:: MMCifInfoStructRefSeqDif

  A particular difference between the deposited sequence and the sequence in 
  the database.

  .. attribute:: rnum

    The residue number (1-based) of the residue in the deposited sequence

    :type: :class:`int`

  .. attribute:: details

    A textual description of the difference, e.g. point mutation, 
    expression tag, purification artifact.

    :type: :class:`str`

.. class:: MMCifInfoRevisions

  Revision history of a PDB entry. If you find a '?' somewhere, this means
  'not set'.

  .. attribute:: date_original

    The date when this entry was seen in PDB for the very first time. This is
    not necessarily the release date. Expected format 'yyyy-mm-dd'.

    :type: :class:`str`

  .. attribute:: first_release

    Index + 1 of the revision releasing this entry. If the value is 0, was not
    set yet. Set first time we encounter a :meth:`GetStatus` value of
    "full release" (mmCIF versions < 5) or "Initial release" (current mmCIF).

    :type: :class:`int`

  .. method:: AddRevision(num, date, status, major=-1, minor=-1)

    Add a new iteration to the history.

    :param num: See :meth:`GetNum`
    :type num:  :class:`int`
    :param date: See :meth:`GetDate`
    :type date:  :class:`str`
    :param status: See :meth:`GetStatus`
    :type status:  :class:`str`
    :param major: See :meth:`GetMajor`
    :type major:  :class:`int`
    :param minor: See :meth:`GetMinor`
    :type minor:  :class:`int`

    :raises: Exception if *num* is <= the last added iteration.

  .. method:: GetSize()

    :return: Number of revisions (valid revision indices are in [0, number-1]).
    :rtype:  :class:`int`

  .. method:: GetDate(i)

    :param i: Index of revision
    :type i:  :class:`int`
    :return: Date the PDB revision took place. Expected format 'yyyy-mm-dd'.
    :rtype:  :class:`str`
    :raises: Exception if *i* out of bounds.

  .. method:: GetNum(i)

    :param i: Index of revision
    :type i:  :class:`int`
    :return: Unique identifier of revision (assigned in increasing order)
    :rtype:  :class:`int`
    :raises: Exception if *i* out of bounds.

  .. method:: GetStatus(i)

    :param i: Index of revision
    :type i: :class:`int`
    :return: The status of this revision.
    :rtype:  :class:`str`
    :raises: Exception if *i* out of bounds.

  .. method:: GetMajor(i)

    :param i: Index of revision
    :type i: :class:`int`
    :return: The major version of this revision (-1 if not set).
    :rtype:  :class:`int`
    :raises: Exception if *i* out of bounds.

  .. method:: GetMinor(i)

    :param i: Index of revision
    :type i: :class:`int`
    :return: The minor version of this revision (-1 if not set).
    :rtype:  :class:`int`
    :raises: Exception if *i* out of bounds.

  .. method:: GetLastDate()

    :return: Date of the latest revision ('?' if no revision set).
    :rtype:  :class:`str`

  .. method:: GetLastMajor()

    :return: Major version of the latest revision (-1 if not set).
    :rtype:  :class:`int`

  .. method:: GetLastMinor()

    :return: Minor version of the latest revision (-1 if not set).
    :rtype:  :class:`int`

  .. method:: SetDateOriginal(date)
              GetDateOriginal()

    See :attr:`date_original`

  .. method:: GetFirstRelease()

    See :attr:`first_release`

.. class:: MMCifInfoEntityBranchLink

  Data from ``pdbx_entity_branch``, most specifically
  ``pdbx_entity_branch_link``. That is connectivity information for branched
  entities, e.g. carbohydrates/ oligosaccharides.
  :class:`Conop Processors <ost.conop.Processor>` can not easily connect them so
  we use this information in :meth:`LoadMMCIF` to do that.

  .. attribute:: atom1

    The first atom of the bond. Corresponds to ``entity_branch_link.atom_id_1``,
    ``entity_branch_link.comp_id_1`` and
    ``entity_branch_link.entity_branch_list_num_1``. Also available via
    :meth:`GetAtom1` and :meth:`SetAtom1`.

    :type: :class:`~ost.mol.AtomHandle`

  .. attribute:: atom2

    The second atom of the bond. Corresponds to ``entity_branch_link.atom_id_2``,
    ``entity_branch_link.comp_id_2`` and
    ``entity_branch_link.entity_branch_list_num_2``. Also available via
    :meth:`GetAtom2` and :meth:`SetAtom2`.

    :type: :class:`~ost.mol.AtomHandle`

  .. attribute:: bond_order

    Order of a bond (e.g. 1=single, 2=double, 3=triple). Corresponds to
    ``entity_branch_link.value_order``. Also available via :meth:`GetBondOrder`
    and :meth:`SetBondOrder`.

    :type: :class:`int`

  .. method:: ConnectBranchLink(editor)

    Establish a bond between :attr:`atom1` and :attr:`atom2` of a
    :class:`MMCifInfoEntityBranchLink`.

    :param editor: The editor instance to call for connecting the atoms.
    :type editor: :class:`~ost.mol.XCSEditor`
    :returns: Nothing

  .. method:: GetAtom1

    See :attr:`atom1`

  .. method:: GetAtom2

    See :attr:`atom2`

  .. method:: GetBondOrder

    See :attr:`bond_order`

  .. method:: SetAtom1

    See :attr:`atom1`

  .. method:: SetAtom2

    See :attr:`atom2`

  .. method:: SetBondOrder

    See :attr:`bond_order`

..  LocalWords:  cas isbn pubmed asu seqres conop ConnectAll casp COMPND OBSLTE
..  LocalWords:  SPRSDE pdb func autofunction exptl attr pdbx oper conf spr dif
..  LocalWords:  biounits biounit uniprot UNP seqs AddMMCifPDBChainTr cif asym
..  LocalWords:  auth GetMMCifPDBChainTr AddPDBCMMCifhainTr GetPDBMMCifChainTr
..  LocalWords:  GetRevisions AddRevision SetRevisionsDateOriginal GetSize str
..  LocalWords:  GetNum num GetStatus GetLastDate GetFirstRelease storable
..  LocalWords:  cas isbn pubmed asu seqres conop casp COMPND OBSLTE LoadMMCIF
..  LocalWords:  SetChainList MMCifInfoTransOp ChainTypes MMCifInfoStructRef
..  LocalWords:  MMCifInfoRevisions bool difs MMCifInfoStructRefSeqDif rnum
..  LocalWords:  SetDateOriginal GetDateOriginal yyyy operationsintervalls
..  LocalWords:  chainintervalls GetChainIntervalList GetMethodDetails GetAtom
..  LocalWords:  GetOperationsIntervalList SetMethodDetails oligosaccharides
..  LocalWords:  SetAtom GetBondOrder SetBondOrder MMCifInfoEntityBranchLink
..  LocalWords:  GetEntityBranchByChain param
