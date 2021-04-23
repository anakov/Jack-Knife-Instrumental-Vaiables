/*
**  setseed(iSeed)
**
**  Purpose:
**    Reset the iSeed, based on iSeed, or on the current time if
**    iSeed <= 0
**
**  Inputs:
**    iSeed       Integer, value for ranseed if iSeed > 0, else current
**                time is used
**
**  Outputs:
**    none        ranseed is reset
**
**  Author:
**    Charles Bos
**
**  Date:
**    24/7/2000
*/
setseed(iSeed)
{
decl h, m, s, stim;

stim = time();

if (iSeed <= 0)
  {
    sscan(stim[0:1], "%d", &h);
    sscan(stim[3:4], "%d", &m);
    sscan(stim[6:7], "%d", &s);
    iSeed= h * 3600 + m * 60 + s;
  }
ranseed(iSeed);
}

SetSeed(iSeed)
{
  setseed(iSeed);
}
