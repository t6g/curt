
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

cxfrt: cxfrt.o  chkopts
	-${CLINKER} -o cxfrt cxfrt.o  ${PETSC_KSP_LIB}
	${RM} cxfrt.o

