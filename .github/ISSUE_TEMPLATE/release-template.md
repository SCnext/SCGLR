---
name: Release template
about: Template for releasing a new version
title: ''
labels: ''
assignees: billy34

---

Prepare for release:

- [ ] git pull
- [ ] Check current CRAN check results
- [ ] Polish NEWS
- [ ] devtools::build_readme()
- [ ] urlchecker::url_check()
- [ ] devtools::check(remote = TRUE, manual = TRUE)
- [ ] devtools::check_win_devel()
- [ ] rhub::check_for_cran()
- [ ] rhub::check(platform = 'ubuntu-rchk')
- [ ] rhub::check_with_sanitizers()
- [ ] revdepcheck::cloud_check()
- [ ] Update cran-comments.md
- [ ] git push

Submit to CRAN:

- [ ] usethis::use_version('patch')
- [ ] devtools::submit_cran()
- [ ] Approve email

Wait for CRAN...

- [ ] Accepted 🎉
- [ ] git push
- [ ] usethis::use_github_release()
- [ ] usethis::use_dev_version()
- [ ] git push
