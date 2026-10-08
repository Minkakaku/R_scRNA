# 按项目清单整理已有结果；正式图仅复制，不重算；原目录完整移入可恢复归档。
function Invoke-ProjectResultsOrganization {
    param([Parameter(Mandatory)][string]$ProjectDir, [string]$PlanFile = 'config/results_organization.json')
    $ErrorActionPreference = 'Stop'
    $projectRoot = (Resolve-Path -LiteralPath $ProjectDir).Path.TrimEnd('\')
    function ProjectPath([string]$relative) {
        $full = [IO.Path]::GetFullPath((Join-Path $projectRoot $relative))
        if (-not $full.StartsWith($projectRoot + '\', [StringComparison]::OrdinalIgnoreCase)) {
            throw "Path outside project: $full"
        }
        return $full
    }
    $plan = Get-Content -LiteralPath (ProjectPath $PlanFile) -Raw | ConvertFrom-Json
    $archive = ProjectPath $plan.archive_directory
    $cache = ProjectPath $plan.cache_directory
    if ((Test-Path -LiteralPath $archive) -or (Test-Path -LiteralPath $cache)) { throw 'Archive or cache already exists' }
    $records = @()
    # 先校验全部源、目的路径和哈希，禁止覆盖已有正式文件。
    foreach ($directory in $plan.source_directories) {
        $source = ProjectPath "results/$directory"
        if (-not (Test-Path -LiteralPath $source -PathType Container)) { throw "Missing source: $source" }
        foreach ($file in Get-ChildItem -LiteralPath $source -Recurse -File) {
            $relative = [IO.Path]::GetRelativePath((ProjectPath 'results'), $file.FullName).Replace('\', '/')
            $records += [pscustomobject]@{ relative = $relative; sha256 = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash; extension = $file.Extension }
        }
    }
    foreach ($item in $plan.exports) {
        $source = ProjectPath "results/$($item.source)"
        $target = ProjectPath "results/$($item.destination)"
        if (-not (Test-Path -LiteralPath $source -PathType Leaf)) { throw "Missing export: $source" }
        if (Test-Path -LiteralPath $target) { throw "Refusing to overwrite: $target" }
    }
    if (@($plan.exports.destination | Select-Object -Unique).Count -ne @($plan.exports).Count) { throw 'Duplicate export destinations' }
    New-Item -ItemType Directory -Path $archive, $cache | Out-Null
    Copy-Item -LiteralPath (ProjectPath 'config/config.yaml') -Destination (Join-Path $archive 'config.before.yaml')
    foreach ($item in $plan.exports) {
        $target = ProjectPath "results/$($item.destination)"
        New-Item -ItemType Directory -Force -Path (Split-Path $target) | Out-Null
        Copy-Item -LiteralPath (ProjectPath "results/$($item.source)") -Destination $target
    }
    # 表达证据、参数、清单和绘图数据进入计算记录目录；结果页不堆积过程表。
    foreach ($record in $records | Where-Object { $_.extension -notin '.png', '.pdf', '.md' }) {
        $target = ProjectPath "$($plan.cache_directory)/$($record.relative)"
        New-Item -ItemType Directory -Force -Path (Split-Path $target) | Out-Null
        Copy-Item -LiteralPath (ProjectPath "results/$($record.relative)") -Destination $target
        if ((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash -ne $record.sha256) { throw "Cache hash mismatch: $target" }
    }
    # 目录移动前明确核验源和目的绝对路径均位于同一项目；原结果字节不改写。
    foreach ($directory in $plan.source_directories) {
        $source = ProjectPath "results/$directory"
        $target = ProjectPath "$($plan.archive_directory)/results/$directory"
        New-Item -ItemType Directory -Force -Path (Split-Path $target) | Out-Null
        Move-Item -LiteralPath $source -Destination $target
    }
    foreach ($record in $records) {
        $target = ProjectPath "$($plan.archive_directory)/results/$($record.relative)"
        if ((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash -ne $record.sha256) { throw "Archive hash mismatch: $target" }
    }
    foreach ($item in $plan.exports) {
        $expected = ($records | Where-Object relative -eq $item.source).sha256
        $target = ProjectPath "results/$($item.destination)"
        if ((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash -ne $expected) { throw "Export hash mismatch: $target" }
    }
    $projectReadme = ProjectPath 'README.md'
    Move-Item -LiteralPath $projectReadme -Destination (ProjectPath "$($plan.archive_directory)/project_README.md")
    foreach ($directory in $plan.remove_empty_directories) {
        $target = ProjectPath "results/$directory"
        if ((Test-Path -LiteralPath $target) -and @(Get-ChildItem -LiteralPath $target -Force).Count -eq 0) {
            Remove-Item -LiteralPath $target
        }
    }
    $configFile = ProjectPath 'config/config.yaml'
    $configText = [IO.File]::ReadAllText($configFile)
    foreach ($directory in $plan.source_directories) {
        $configText = $configText.Replace("results/$directory", "$($plan.cache_directory)/$directory")
    }
    [IO.File]::WriteAllText($configFile, $configText, [Text.UTF8Encoding]::new($false))
    $records | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath (ProjectPath $plan.integrity_log) -Encoding utf8
    Write-Output "RESULTS_ORGANIZATION_OK archived=$($records.Count) published=$(@($plan.exports).Count) hashes=all_matched"
}
