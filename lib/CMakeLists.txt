file(GLOB COND_FRMW_SRC *.cc enumeration_algorithm/*.cc)

add_library(conditions_algorithms
        SHARED ${COND_FRMW_SRC})
