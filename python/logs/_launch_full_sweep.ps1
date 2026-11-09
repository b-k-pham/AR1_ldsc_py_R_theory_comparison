$env:CONDA_PREFIX = 'C:\Users\Ben Pham\.conda\envs\gwash'
$env:PATH = 'C:\Users\Ben Pham\.conda\envs\gwash;C:\Users\Ben Pham\.conda\envs\gwash\Scripts;C:\Users\Ben Pham\.conda\envs\gwash\Library\bin;' + $env:PATH
Set-Location 'D:\PhD\gwash_sim_test\to_github_restore\compare_implementations\python'
.\run_n10000_other3.ps1 -ThreadLimit 10
