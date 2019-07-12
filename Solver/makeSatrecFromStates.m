% make Satrec from states
function satrec = makeSatrecFromStates(satrec,Xest)
satrec.ecco = Xest(1);
satrec.inclo = Xest(2);
satrec.nodeo = Xest(3);
satrec.argpo = Xest(4);
satrec.mo = Xest(5);
satrec.no_kozai = Xest(6);
satrec.bstar = Xest(7);
end