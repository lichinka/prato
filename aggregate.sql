--UPDATE coverage_CGBRIA SET rscp=27.2-rscp;
--UPDATE coverage_CGBRIB SET rscp=27.2-rscp;
--UPDATE coverage_CGBRIC SET rscp=27.2-rscp;
--UPDATE coverage_CSOSTA SET rscp=27.0-rscp;
--UPDATE coverage_CSOSTB SET rscp=28.0-rscp;
--UPDATE coverage_CSOSTC SET rscp=30.0-rscp;

SELECT east, north, MAX(rscp)
  FROM (SELECT east, north, rscp FROM coverage_CGBRIA
        UNION
        SELECT east, north, rscp FROM coverage_CGBRIB
        UNION
        SELECT east, north, rscp FROM coverage_CGBRIC
        UNION
        SELECT east, north, rscp FROM coverage_CSOSTA
        UNION
        SELECT east, north, rscp FROM coverage_CSOSTB
        UNION
        SELECT east, north, rscp FROM coverage_CSOSTC) partial
GROUP BY partial.east, partial.north;
