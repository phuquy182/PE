# Kiểm nhánh lắp ghép bằng phần tử hữu hạn (việc 3, Hướng 1, Bài báo 2)

*Ngày 13/11/2026. Ngoài đề cương luận văn: chỉ dùng cho Bài báo 2, **chờ thầy duyệt**. Không đưa vào luận văn.*

Chuyển giao: `claude/ChuyenGiao_KiemFE_Huong1_20261113.md`. Báo cáo PDF: `BaoCao_KiemFE_Huong1_20261113.pdf`. Mã FE: thư mục `fe_fitcheck/` (scikit-fem 12.0.2 + triangle).

## 1. Kết luận

- Ca A (vòng trơn D_eff theo [9]) trùng nghiệm Lamé: sai lệch lớn nhất 0,0009 % → mã đúng.
- Ca C, D_eff theo [9]: σθ,bore −0,8 … −0,2 %; σθ,z0 −1,7 … −1,5 %; k_cl = d/D_eff −1,6 … −1,5 % → **khớp** (≤ 5 %).
- Ca C, p_f trung bình: công thức [9] thấp hơn FE 10,4 … 13,6 % → **không khớp** (> 10 %). Tỉ số FE/công thức 1,117 … 1,157.
- D_eff cùng diện tích: p_f −1,7 … −1,6 % (khớp); σθ,bore −3,2 … −2,9 %, σθ,z0 −5,0 … −4,0 %, k_cl −5,3 … −4,2 %. Không có D_eff nào khớp cả bốn đại lượng.
- Chiều sai lệch với tuổi thọ: σθ,z0 công thức thấp hơn FE → hơi thiên về không an toàn. Ca 6310 (5 200 N, 18 000 h, 1 500 v/ph): nhân K_δ với 1,0173 → δ_life 15,70 → 15,50 µm. (Chỉ minh hoạ, không sửa chương trình.)
- p_f không đi thẳng vào chuỗi tuổi thọ (σθ = p_f·g(D_eff); δ_slip theo NTN [4]). Sai lệch p_f chỉ ảnh hưởng tới giá trị p_f được in ra.
- Thăm dò vai ±0,05·D_w: σθ,z0 đổi ≤ 0,23 %, k_cl ≤ 0,29 %, p_f ≤ 4,6 %.

## 2. Định nghĩa D_eff trong [9] (đã mở toàn văn)

[9] = NASA/TM-2010-216913. Tr. 4, pt. (6) – (8a): θ = 1,854 rad, A_sector = (R²/2)(θ − sin θ) = 0,1118·d² (d: đường kính bi), D_eff = D_IR + 0,2·d − A_sector/W. Văn bản ghi: D_eff tính «bằng cách cộng diện tích vai vào vòng không vai rồi chia cho bề rộng». Hình 6, tr. 5 ghi chiều cao vai 0,2·d theo bán kính. Ví dụ ổ cỡ 210, bước (2), tr. 16: D_eff = 58,94 mm (trùng `B6210.D_eff`). Pt. (8a) cộng (0,2·d − A/W) một lần vào đường kính, tức bằng nửa phần tăng bán kính của vòng trơn cùng diện tích (D_IR + 0,4·d − 2A/W). Suy luận ngày 13/11 khớp nguyên văn.

## 3. Bảng chính: ca C, δ_eff = 10 µm (sai lệch = (công thức − FE)/FE)

| Ổ | D_eff [9] | D_eff cùng d.t. | p_f FE | p_f [9] | p_f c.d.t. | σθ,bore FE | [9] | c.d.t. | σθ,z0 FE | [9] | c.d.t. | k_cl FE | [9] | c.d.t. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 6306 | 39,89 | 41,48 | 16,65 | 14,91 (−10,4 %) | 16,37 (−1,7 %) | 53,84 | 53,72 (−0,2 %) | 52,26 (−2,9 %) | 41,40 | 40,77 (−1,5 %) | 39,31 (−5,0 %) | 0,764 | 0,752 (−1,5 %) | 0,723 (−5,3 %) |
| 6210 | 58,94 | 60,58 | 6,68 | 5,77 (−13,6 %) | 6,56 (−1,7 %) | 35,70 | 35,40 (−0,8 %) | 34,61 (−3,0 %) | 31,11 | 30,66 (−1,5 %) | 29,86 (−4,0 %) | 0,861 | 0,848 (−1,5 %) | 0,825 (−4,2 %) |
| 6310 | 63,26 | 65,57 | 8,76 | 7,73 (−11,8 %) | 8,61 (−1,6 %) | 33,62 | 33,45 (−0,5 %) | 32,56 (−3,2 %) | 27,37 | 26,91 (−1,7 %) | 26,02 (−4,9 %) | 0,803 | 0,790 (−1,6 %) | 0,763 (−5,1 %) |

Đơn vị: mm, MPa. σθ,bore tại mặt phẳng giữa. k_cl FE = 2·u_r(đáy rãnh, z = 0)/δ_eff; công thức = d/D_eff.

D tính ngược từ p_f FE (D làm pt. (3.9) cho đúng p_f FE): 6306: 41,81 mm; 6210: 60,83 mm; 6310: 65,96 mm.

## 4. Kiểm mô hình

- Hội tụ lưới (6310 ca C, ba mức, cỡ phần tử dưới rãnh z0/2, z0/4, z0/8 với z0 = 0,2077 mm): σθ,z0 = 27,3721; 27,3724; 27,3725 MPa → đổi 0,0004 % giữa hai mức cuối (tiêu chí ≤ 1 %).
- Tuyến tính: δ 30/10 µm → tỉ số p_f 3,00000, σθ,z0 3,00000.
- Trục dài thêm 2B thay vì B: p_f +0,04 %, σθ,z0 +0,02 %.
- Cân bằng: ∫σθ dA trên mặt cắt nửa vòng = p_f·(d/2)(B/2), sai lệch < 1e-11 (6306, 6310, ca B và C).
- Tiếp xúc: mọi nhân tử Lagrange đều nén ở mọi ca.
- Gói `bearing_fit` dựng lại từ Drive (27/10) đạt hồi quy: f_fit = 0,8591; δ_max,life = 15,70 µm.

## 5. Ca B và D (sai lệch của công thức [9] so với FE)

| Ổ | Ca | vai·D_w | p_f FE | Δp_f | σθ,z0 FE | Δσθ,z0 | k_cl FE | Δk_cl |
|---|---|---|---|---|---|---|---|---|
| 6306 | B | 0,20 | 16,217 | −8,1 % | 42,68 | −4,5 % | 0,7661 | −1,8 % |
| 6306 | C | 0,20 | 16,646 | −10,4 % | 41,40 | −1,5 % | 0,7635 | −1,5 % |
| 6306 | D | 0,15 | 16,109 | −7,5 % | 41,32 | −1,3 % | 0,7655 | −1,8 % |
| 6306 | D | 0,25 | 17,090 | −12,8 % | 41,49 | −1,7 % | 0,7622 | −1,3 % |
| 6210 | B | 0,20 | 6,506 | −11,3 % | 31,43 | −2,5 % | 0,8537 | −0,6 % |
| 6210 | C | 0,20 | 6,678 | −13,6 % | 31,11 | −1,5 % | 0,8613 | −1,5 % |
| 6210 | D | 0,15 | 6,371 | −9,4 % | 31,11 | −1,5 % | 0,8637 | −1,8 % |
| 6210 | D | 0,25 | 6,942 | −16,9 % | 31,12 | −1,5 % | 0,8593 | −1,3 % |
| 6310 | B | 0,20 | 8,524 | −9,4 % | 28,01 | −3,9 % | 0,8025 | −1,5 % |
| 6310 | C | 0,20 | 8,756 | −11,8 % | 27,37 | −1,7 % | 0,8034 | −1,6 % |
| 6310 | D | 0,15 | 8,437 | −8,4 % | 27,33 | −1,6 % | 0,8054 | −1,9 % |
| 6310 | D | 0,25 | 9,017 | −14,3 % | 27,42 | −1,9 % | 0,8020 | −1,4 % |

## 6. Giới hạn

- Vai 0,2·D_w và biên dạng vai là giả định; bỏ mép vát. Không ma sát, trục đặc, đàn hồi tuyến tính.
- Đỉnh áp suất ở mép bề mặt lắp (ca C) là kỳ dị, phụ thuộc lưới; chỉ dùng p_f trung bình.
- FE không kiểm chứng tuổi thọ, không thay thực nghiệm.

## 7. Việc chờ

1. Thầy duyệt việc đưa mục kiểm FE vào Bài báo 2.
2. p_f theo [9] lệch FE trên 10 %: Bài báo 1 và luận văn (ví dụ Bảng 4.3) có in giá trị p_f không? Tác giả quyết định cách xử lý; chưa sửa gì.
3. Tác giả có bản vẽ hoặc số đo vai, mép vát của 6306, 6210, 6310 không? Có thì chạy lại ca C.
4. Chép `fe_fitcheck/` vào `07_BaiBao\BaiBao2_QuocTe\` (thư mục con mới) khi máy kết nối, báo trước khi chép.

## 8. Tệp

- `fe_fitcheck/fe_axisym.py`: lưới, lời giải, hậu xử lý.
- `fe_fitcheck/run_fe_check.py`: mọi ca → `results/fe_results.json`, `.csv`.
- `fe_fitcheck/make_figures.py`, `make_report.py`, `make_note.py`: hình, PDF, ghi chú này.
- `bearing_fit/`: bản dựng lại để chạy kiểm (không sửa nội dung).
