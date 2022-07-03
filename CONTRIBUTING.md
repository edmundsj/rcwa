# How to contribute
First of all, thank you for taking the time to contribute to this project! I've tried to make a stable project and to fix bugs and add features as I am able. :)

## Getting started

### Check out the roadmap

I have some functionalities in mind and we have issued them and there is a *milestone* label available on the issue. If there is a bug or a feature that is not listed in the **issues** page or there is no one assigned to the issue, feel free to fix/add it! Although it's better to discuss it in the issue or create a new issue for it so there is no confilcting code.

### Writing some code!

Contributing to a project on Github is pretty straight forward. If this is you're first time, these are the steps you should take.

- Fork this repo.

And that's it! Read the code available and change the part you don't like! You're change should not break the existing code and should pass the tests.

If you're adding a new functionality, start from the branch **master**. It would be a better practice to create a new branch and work in there. The **develop** branch contains features in-progress which are locally stable, so you may want to check that out to see what is being worked on. Other feature-specific branches may exist as well.

When you're done, submit a pull request and I will check it out. I'll let you know if there is any problem or any changes that should be considered. My main comments will probably have to do with testing and documentation (see below).

### Tests

I've written a comprehensive test suite and you can run it to assure the stability of the code, just try `pytest` in the package folder. This runs the same test suite that runs on the build server. This package aims to be production-ready, so any new features *must* be thoroughly tested before being merged. If you are working on a new feature (first, awesome :D) please write unit and integration tests for it.

### Documentation

Every chunk of code that may be hard to understand has some comments above it. If you write some new code or change some part of the existing code in a way that it would not be functional without changing it's usages, it needs to be documented.
