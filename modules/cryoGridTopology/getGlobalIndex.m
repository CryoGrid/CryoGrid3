function idx=getGlobalIndex( mesoIndex, microIndex, Nmicro )
    assert( mesoIndex<=length(Nmicro), 'mesoIndex out of range' );
    assert( microIndex<=Nmicro(mesoIndex), 'micorIndex out of range');
    idx=sum(Nmicro(1:mesoIndex-1))+microIndex;
end