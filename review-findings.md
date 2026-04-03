## REVIEW CLEAN
## Code Review: sheaf_nma.html
### Date: 2026-04-03
### Summary: 1 P0 FIXED, 2 P1, 3 P2

#### P0 -- Critical (FIXED)
- **[FIXED] P0-1** [Security]: CSV export (`exportCSV`) writes study names directly into CSV without formula injection protection. Custom CSV input allows arbitrary study names starting with `=`, `+`, `@`, `\t`, `\r` -- these can trigger formula execution when opened in Excel/Sheets. Fix: add `sanitizeCSVCell()` to prepend `'` to cells starting with dangerous characters (per rule: NOT `-`).

#### P1 -- Important
- **P1-1** [Accessibility]: Test panel (`#testPanel`) has no Escape key handler. The TDA_MA sibling project correctly implements `keydown`/`removeEventListener` for its test modal, but SheafNMA's bottom panel relies only on a close button click. Low risk since the panel is not a modal overlay.
- **P1-2** [Accessibility]: Tab buttons use `data-tab="data"` but the handler constructs panel ID as `"tab-" + "data"` = `"tab-data"`. This works correctly since panel IDs match (`tab-data`, `tab-network`, etc.), but the convention is inconsistent with typical `aria-controls` patterns.

#### P2 -- Minor
- **P2-1** [Stats]: Sheaf Laplacian coboundary convention uses `effect = x[treat2] - x[treat1]` which is consistent with the graph-theoretic orientation. The node estimation via reduced system (fixing node 0 = 0) is standard NMA practice. No issues found.
- **P2-2** [Stats]: GII formula `totalResidSS / df` where `df = m - (n-1)` is the correct Q-statistic analogue for network consistency. Matches Higgins/Jackson design-by-treatment interaction test.
- **P2-3** [Code Quality]: `escapeHtml` correctly escapes all 5 characters (`<`, `>`, `&`, `"`, `'`), verified in test 14. All innerHTML paths use `escapeHtml()` for user data.

#### Structural Checks
- Div balance: 64 open, 64 close -- balanced
- `</script>` literal: only at line 2548 (actual closing tag) -- safe
- `</html>` closing tag: present at line 2550
- PRNG: uses `xoshiro128**` with seed -- deterministic
- Blob URLs: `URL.revokeObjectURL()` called after download -- no leak

#### Test Results: 24/24 pass (in-browser suite)
