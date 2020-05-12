function [meso,micro]=getMesoMicroIndices( globalIndex, Nmicro )
    assert( globalIndex<=sum(Nmicro), 'globalIndex out of range');
    cs=cumsum( Nmicro );
    meso=find( cs >= globalIndex, 1, 'first');
    micro=globalIndex-sum(Nmicro(1:meso-1));
end