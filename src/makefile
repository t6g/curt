
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

cxrt: cxrt.o  chkopts
	-${CLINKER} -o cxrt cxrt.o  ${PETSC_KSP_LIB}
	${RM} cxrt.o

