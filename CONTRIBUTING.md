# How to Contribute
## Step 1:
After cloning the repo, create a new branch with:
```
git checkout -b branch_name
```
Name the branch something reflective of what you're working on.

## Step 2:
Program! Implement the feature, make sure it works free of bugs and memory leaks. Document the the code with at least function banners.
Feel free to commit and push your code (working or not), to your own branch.
When you push your work to the repo for the first time, you have to tell github that you have made a new branch. Push using:
```
git push --set-upstream origin branch_name
```
Once you have made this first push, you no longer need to use `--set-upstream origin` when you push going forward.

## Step 3:
Merge your work to the main branch. To do this in the safest way possible, start by executing:
```
git merge main
```
This will update your branch with anything that might be new in main since you've started working on your branch.  

## Step 4a:
If things merge with no issues:
1. Proceed by moving to the main
  - `git checkout main`
2. Merge your finished branch with main
  - `git merge branch_name`
3. Push!
  - `git push`

And you're done!
  
## Step 4b:
If you get merge conflicts:
Resolve the merge conflict. This happens when git isn't sure how to combine the new code and the old code.
It usually looks like this:
```
int main () {

<<<<<<< HEAD
printf("Hello, I was written in the main branch!");
=======
printf("Howdy, I am the new feature from the branch_name branch");
>>>>>>> branch_name

  return 0;
}
```
It might not always be as simple as choosing between the top and bottom, sometimes you must combine code from each section.
Combine the two sections appropiately, recomile, and make sure everything is working as expected. Fix it otherwise.
Once working, commit with a message telling us about the merge, then continue at step 4a.

If the merge conflict gets too out of hand, you can bail on the merge with:
```
git reset --hard HEAD
```
Then reach out to @m18coppola for help.
