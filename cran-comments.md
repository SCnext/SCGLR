- Ill formed urls not caught by below test were fixed.
- Fixed author differing from Author@R (documentation states that comment field 
of person() is free but in fact not! Only two names are allowed ORCID and ROR 
other string must not have names)

This is an update of the package introducing new functionalities 
and providing bug fixes.

## Test environments
* win-builder
* rhub windows-latest
* rhub ubuntu-next

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

No package depends on this one and only additions so should be backward compatible.
