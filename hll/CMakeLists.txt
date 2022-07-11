set(HLL_SOURCES hll.c hll_empirical.c)
set_source_files_properties(${HLL_SOURCES} PROPERTIES LANGUAGE C)

add_library(hll STATIC ${HLL_SOURCES})
target_compile_options(hll PRIVATE -Wall -Wextra -Werror -Wpedantic)
target_include_directories(hll INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")
