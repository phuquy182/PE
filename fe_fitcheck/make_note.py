# -*- coding: utf-8 -*-
"""make_note.py — Writes the project note claude/KiemFE_Huong1_20261113.md from the results."""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_report as R   # noqa: E402  (reuses the parsed results and derived numbers)

vn, vp, pct, row, FORM, NAMES = R.vn, R.vp, R.pct, R.row, R.FORM, R.NAMES


def main(path):
    L = []
    a = L.append
    a('# Kiểm nhánh lắp ghép bằng phần tử hữu hạn (việc 3, Hướng 1, Bài báo 2)\n')
    a('*Ngày 13/11/2026. Ngoài đề cương luận văn: chỉ dùng cho Bài báo 2, **chờ thầy duyệt**. '
      'Không đưa vào luận văn.*\n')
    a('Chuyển giao: `claude/ChuyenGiao_KiemFE_Huong1_20261113.md`. Báo cáo PDF: '
      '`BaoCao_KiemFE_Huong1_20261113.pdf`. Mã FE: thư mục `fe_fitcheck/` (scikit-fem 12.0.2 + triangle).\n')
    a('## 1. Kết luận\n')
    a(f'- Ca A (vòng trơn D_eff theo [9]) trùng nghiệm Lamé: sai lệch lớn nhất {vn(R.errA, 4)} % → mã đúng.')
    a(f'- Ca C, D_eff theo [9]: σθ,bore {vp(R.sb9[0])} … {vp(R.sb9[1])} %; σθ,z0 {vp(R.sz9[0])} … '
      f'{vp(R.sz9[1])} %; k_cl = d/D_eff {vp(R.k9[0])} … {vp(R.k9[1])} % → **khớp** (≤ 5 %).')
    a(f'- Ca C, p_f trung bình: công thức [9] thấp hơn FE {vn(-R.p9[1], 1)} … {vn(-R.p9[0], 1)} % → '
      f'**không khớp** (> 10 %). Tỉ số FE/công thức {vn(min(R.ratio_p), 3)} … {vn(max(R.ratio_p), 3)}.')
    a(f'- D_eff cùng diện tích: p_f {vp(R.pea[0])} … {vp(R.pea[1])} % (khớp); σθ,bore {vp(R.sbea[0])} … '
      f'{vp(R.sbea[1])} %, σθ,z0 {vp(R.szea[0])} … {vp(R.szea[1])} %, k_cl {vp(R.kea[0])} … {vp(R.kea[1])} %. '
      'Không có D_eff nào khớp cả bốn đại lượng.')
    a(f'- Chiều sai lệch với tuổi thọ: σθ,z0 công thức thấp hơn FE → hơi thiên về không an toàn. Ca 6310 '
      f'(5 200 N, 18 000 h, 1 500 v/ph): nhân K_δ với {vn(R.ratio_z["6310"], 4)} → δ_life {vn(R.dl0)} → '
      f'{vn(R.dl1)} µm. (Chỉ minh hoạ, không sửa chương trình.)')
    a('- p_f không đi thẳng vào chuỗi tuổi thọ (σθ = p_f·g(D_eff); δ_slip theo NTN [4]). Sai lệch p_f chỉ '
      'ảnh hưởng tới giá trị p_f được in ra.')
    a(f'- Thăm dò vai ±0,05·D_w: σθ,z0 đổi ≤ {vn(R.shoulder_dz, 2)} %, k_cl ≤ {vn(R.shoulder_dk, 2)} %, '
      f'p_f ≤ {vn(R.shoulder_dp, 1)} %.\n')
    a('## 2. Định nghĩa D_eff trong [9] (đã mở toàn văn)\n')
    a('[9] = NASA/TM-2010-216913. Tr. 4, pt. (6) – (8a): θ = 1,854 rad, A_sector = (R²/2)(θ − sin θ) = '
      '0,1118·d² (d: đường kính bi), D_eff = D_IR + 0,2·d − A_sector/W. Văn bản ghi: D_eff tính «bằng cách cộng '
      'diện tích vai vào vòng không vai rồi chia cho bề rộng». Hình 6, tr. 5 ghi chiều cao vai 0,2·d theo bán '
      'kính. Ví dụ ổ cỡ 210, bước (2), tr. 16: D_eff = 58,94 mm (trùng `B6210.D_eff`). Pt. (8a) cộng '
      '(0,2·d − A/W) một lần vào đường kính, tức bằng nửa phần tăng bán kính của vòng trơn cùng diện tích '
      '(D_IR + 0,4·d − 2A/W). Suy luận ngày 13/11 khớp nguyên văn.\n')
    a('## 3. Bảng chính: ca C, δ_eff = 10 µm (sai lệch = (công thức − FE)/FE)\n')
    a('| Ổ | D_eff [9] | D_eff cùng d.t. | p_f FE | p_f [9] | p_f c.d.t. | σθ,bore FE | [9] | c.d.t. | '
      'σθ,z0 FE | [9] | c.d.t. | k_cl FE | [9] | c.d.t. |')
    a('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for n in NAMES:
        r, f9, fa = row(n, 'C'), FORM[n]['ref9'], FORM[n]['equal_area']
        cells = [n, vn(f9['D_eff']), vn(fa['D_eff'])]
        for kf, kr, nd in (('p_f', 'p_f', 2), ('s_bore', 's_bore_mid', 2), ('s_z0', 's_z0', 2), ('k_cl', 'k_cl', 3)):
            cells += [vn(r[kr], nd), f'{vn(f9[kf], nd)} ({vp(pct(r[kr], f9[kf]))} %)',
                      f'{vn(fa[kf], nd)} ({vp(pct(r[kr], fa[kf]))} %)']
        a('| ' + ' | '.join(cells) + ' |')
    a('\nĐơn vị: mm, MPa. σθ,bore tại mặt phẳng giữa. k_cl FE = 2·u_r(đáy rãnh, z = 0)/δ_eff; công thức = '
      'd/D_eff.\n')
    a('D tính ngược từ p_f FE (D làm pt. (3.9) cho đúng p_f FE): ' + '; '.join(
        f'{n}: {vn(row(n, "C")["D_eff_backcalc_p"])} mm' for n in NAMES) + '.\n')
    a('## 4. Kiểm mô hình\n')
    c = R.conv
    a(f'- Hội tụ lưới (6310 ca C, ba mức, cỡ phần tử dưới rãnh z0/2, z0/4, z0/8 với z0 = {vn(c[0]["z0"], 4)} mm): '
      f'σθ,z0 = {vn(c[0]["s_z0"], 4)}; {vn(c[1]["s_z0"], 4)}; {vn(c[2]["s_z0"], 4)} MPa → đổi '
      f'{vn(abs(pct(c[1]["s_z0"], c[2]["s_z0"])), 4)} % giữa hai mức cuối (tiêu chí ≤ 1 %).')
    a(f'- Tuyến tính: δ 30/10 µm → tỉ số p_f {vn(R.lin["p_f"] / R.c6310["p_f"], 5)}, σθ,z0 '
      f'{vn(R.lin["s_z0"] / R.c6310["s_z0"], 5)}.')
    a(f'- Trục dài thêm 2B thay vì B: p_f {vp(pct(R.c6310["p_f"], R.ext["p_f"]), 2)} %, σθ,z0 '
      f'{vp(pct(R.c6310["s_z0"], R.ext["s_z0"]), 2)} %.')
    a('- Cân bằng: ∫σθ dA trên mặt cắt nửa vòng = p_f·(d/2)(B/2), sai lệch < 1e-11 (6306, 6310, ca B và C).')
    a('- Tiếp xúc: mọi nhân tử Lagrange đều nén ở mọi ca.')
    a('- Gói `bearing_fit` dựng lại từ Drive (27/10) đạt hồi quy: f_fit = 0,8591; δ_max,life = 15,70 µm.\n')
    a('## 5. Ca B và D (sai lệch của công thức [9] so với FE)\n')
    a('| Ổ | Ca | vai·D_w | p_f FE | Δp_f | σθ,z0 FE | Δσθ,z0 | k_cl FE | Δk_cl |')
    a('|---|---|---|---|---|---|---|---|---|')
    for n in NAMES:
        f9 = FORM[n]['ref9']
        for case, sh in (('B', 0.2), ('C', 0.2), ('D', 0.15), ('D', 0.25)):
            r = row(n, case, sh)
            a(f'| {n} | {case} | {vn(sh)} | {vn(r["p_f"], 3)} | {vp(pct(r["p_f"], f9["p_f"]))} % | '
              f'{vn(r["s_z0"])} | {vp(pct(r["s_z0"], f9["s_z0"]))} % | {vn(r["k_cl"], 4)} | '
              f'{vp(pct(r["k_cl"], f9["k_cl"]))} % |')
    a('\n## 6. Giới hạn\n')
    a('- Vai 0,2·D_w và biên dạng vai là giả định; bỏ mép vát. Không ma sát, trục đặc, đàn hồi tuyến tính.')
    a('- Đỉnh áp suất ở mép bề mặt lắp (ca C) là kỳ dị, phụ thuộc lưới; chỉ dùng p_f trung bình.')
    a('- FE không kiểm chứng tuổi thọ, không thay thực nghiệm.\n')
    a('## 7. Việc chờ\n')
    a('1. Thầy duyệt việc đưa mục kiểm FE vào Bài báo 2.')
    a(f'2. p_f theo [9] lệch FE trên 10 %: Bài báo 1 và luận văn (ví dụ Bảng 4.3) có in giá trị p_f không? Tác '
      'giả quyết định cách xử lý; chưa sửa gì.')
    a('3. Tác giả có bản vẽ hoặc số đo vai, mép vát của 6306, 6210, 6310 không? Có thì chạy lại ca C.')
    a('4. Chép `fe_fitcheck/` vào `07_BaiBao\\BaiBao2_QuocTe\\` (thư mục con mới) khi máy kết nối, báo trước khi chép.\n')
    a('## 8. Tệp\n')
    a('- `fe_fitcheck/fe_axisym.py`: lưới, lời giải, hậu xử lý.')
    a('- `fe_fitcheck/run_fe_check.py`: mọi ca → `results/fe_results.json`, `.csv`.')
    a('- `fe_fitcheck/make_figures.py`, `make_report.py`, `make_note.py`: hình, PDF, ghi chú này.')
    a('- `bearing_fit/`: bản dựng lại để chạy kiểm (không sửa nội dung).')
    Path(path).write_text('\n'.join(L) + '\n', encoding='utf-8')
    print(path)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else HERE.parent / 'claude' / 'KiemFE_Huong1_20261113.md')
