var manual =
[
    [ "Core functionality", "core.html", [
      [ "Basic classes", "core.html#basic_classes", [
        [ "Atom", "core.html#atoms", null ],
        [ "System", "core.html#system", null ],
        [ "Frame", "core.html#frame", null ],
        [ "Selection", "core.html#selection", null ]
      ] ],
      [ "Loading molecular system", "core.html#load", [
        [ "Supported file formats", "core.html#formats", null ],
        [ "What is loaded from data files?", "core.html#types_of_info", null ],
        [ "Simple loading", "core.html#simple_load", null ],
        [ "Advanced loading with file handlers", "core.html#advanced_load", null ],
        [ "Loading with callback", "core.html#load_callback", null ],
        [ "Using input filters", "core.html#load_filters", null ]
      ] ],
      [ "Selecting atoms", "core.html#making_selections", [
        [ "Ways of creating selections", "core.html#sel_methods", null ],
        [ "Arguments of selection methods", "core.html#sel_args", null ],
        [ "Modifying selections", "core.html#sel_modify", null ],
        [ "Selecting everything", "core.html#all_sel", null ],
        [ "Selection language", "core.html#sel_lang", [
          [ "Selecting everything", "core.html#select_all", null ],
          [ "Keyword selections", "core.html#keyword_sel", null ],
          [ "Numeric properties", "core.html#num_prop", null ],
          [ "Numeric expressions", "core.html#num_expr", null ],
          [ "Numeric comparisons", "core.html#num_comp", null ],
          [ "Logical expressions", "core.html#log_expr", null ],
          [ "Within selections", "core.html#within_sel", null ],
          [ "By residue selections", "core.html#by_res", null ],
          [ "Distance selections", "core.html#dist_sel", null ]
        ] ],
        [ "Text-based and coordinate-depensent selections", "core.html#text_based_sel", null ],
        [ "Sub-selections", "core.html#subsel", null ],
        [ "Combining selections", "core.html#sel_comb", null ]
      ] ],
      [ "Accessing properties of selected atoms", "core.html#access", [
        [ "Accessor methods in C++", "core.html#cpp_access", null ],
        [ "Accessor methods in Python", "core.html#py_access", null ],
        [ "Indexing operator of Selection object", "core.html#indexing", null ],
        [ "Iterating over selected atoms", "core.html#sel_iter", null ],
        [ "Getting particular property of all selected atoms", "core.html#uniform_prop", null ]
      ] ],
      [ "Building molecular systems", "core.html#sys_build", [
        [ "Adding and deleting atoms", "core.html#atom_add_del", null ],
        [ "Appending and removing systems and selections", "core.html#append_remove", null ],
        [ "Multiplying selections", "core.html#distib", null ],
        [ "Rearranging atoms", "core.html#rearrange", null ]
      ] ],
      [ "Working with periodicity", "core.html#pbc", [
        [ "Getting, setting and modifying periodic box", "core.html#pbc_get_set", null ],
        [ "Wrapping and unwrapping", "core.html#wrap", null ],
        [ "Periodic distances and closest images", "core.html#pbc_measure", null ]
      ] ],
      [ "Fast distance search and spatial grids", "core.html#dist_search", [
        [ "Different forms of distance search", "core.html#dist_search_variants", [
          [ "Searching contacts within single selection", "core.html#search_1_sel", null ],
          [ "Searching contacts between two selections", "core.html#search_2_sel", null ],
          [ "Searching around given selection", "core.html#search_within", null ],
          [ "Repeated searching around multiple targets", "core.html#repeated_search", null ]
        ] ],
        [ "Custom spatial grids", "core.html#custom_grid", null ]
      ] ],
      [ "Evaluating non-bond energies", "core.html#energy", [
        [ "Generating Pteros .pttop files", "core.html#get_top", null ],
        [ "Energy components", "core.html#en_components", null ],
        [ "Computing non-bond energies", "core.html#get_en", [
          [ "Methods of System class", "core.html#en_system", null ],
          [ "Methods of Selection class", "core.html#en_sel", null ]
        ] ]
      ] ],
      [ "Secondary structure of proteins", "core.html#dssp", null ],
      [ "Solvent accessible surface area", "core.html#measure_sasa", null ]
    ] ],
    [ "Analysis of trajectories", "analysis.html", [
      [ "Asynchronous parallel trajectory processing in Pteros", "analysis.html#anal_concepts", null ],
      [ "Trajectory reader and analysis tasks", "analysis.html#processors_consumers", [
        [ "Serial and parallel tasks", "analysis.html#task_types", [
          [ "Selections in parallel tasks", "analysis.html#par_select", null ]
        ] ],
        [ "Processing command-line options", "analysis.html#options", [
          [ "Value suffixes", "analysis.html#suffix", null ]
        ] ],
        [ "Understanding frame metadata", "analysis.html#frame_info", null ],
        [ "\"Magic variables\" in analysis tasks", "analysis.html#magic_vars", null ],
        [ "Removing jumps over periodic boundaries", "analysis.html#jump", null ]
      ] ],
      [ "Analysis plugins", "analysis.html#plugins", null ]
    ] ],
    [ "Pteros tools", "tools.html", null ]
];