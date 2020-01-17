configure_file(config.h.in config.h @ONLY IMMEDIATE)

install(FILES mcb.hpp fvs.hpp dijkstra.hpp mcb_sva_signed.hpp mcb_sva_signed_parallel.hpp signed_dijkstra.hpp spvecgf2.hpp forestindex.hpp mcb_sva_signed_parallel2.hpp mcb_sva_trees.hpp sptrees.hpp util.hpp DESTINATION include/mcb)
