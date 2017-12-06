from cPosValCalculation cimport *
from libc.stdlib cimport malloc, free

cdef struct ChrPosVal:
    cdef char * chrom
    cdef PosVal * posval

cdef class CbedGraph:
    cdef:
        ChrPosVal *data
        int num_of_chrom
        dict chrom

    def __init__ ( self, dict chrom ):
        self.chrom = chrom
        self.num_of_chrom = len( self.chrom.keys() )

    def __cinit__ ( self ):
        cdef int i
        
        self.data = <ChrPosVal *> malloc( self.num_of_chrom * sizeof( ChrPosVal ) )

        for i in range( self.num_of_chrom ):
            self.data[ i ].chrom = <char *> self.chrom[ i ]

    def __dealloc__ ( self ):
        for i in range( self.num_of_chrom ):
            free( self.data[ i ].chrom )
            free( self.data[ i ].posval )
        free( self.data )

    
    
