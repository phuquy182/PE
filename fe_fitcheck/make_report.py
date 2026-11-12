# -*- coding: utf-8 -*-
"""
make_report.py — Builds the Vietnamese PDF report of the FE check from
results/fe_results.json and the figures. Every number in the report is read from
the result file; nothing is typed by hand.
"""
import json
import sys
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
                                Image, KeepTogether)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from bearing_fit import B6310, delta_max_life   # noqa: E402

RES = HERE / 'results'
import matplotlib
FONT_DIR = Path(matplotlib.get_data_path()) / 'fonts' / 'ttf'
pdfmetrics.registerFont(TTFont('DV', str(FONT_DIR / 'DejaVuSans.ttf')))
pdfmetrics.registerFont(TTFont('DVB', str(FONT_DIR / 'DejaVuSans-Bold.ttf')))
pdfmetrics.registerFont(TTFont('DVI', str(FONT_DIR / 'DejaVuSans-Oblique.ttf')))
pdfmetrics.registerFontFamily('DV', normal='DV', bold='DVB', italic='DVI', boldItalic='DVB')

S_BODY = ParagraphStyle('b', fontName='DV', fontSize=9, leading=12.2, alignment=TA_LEFT,
                        spaceAfter=3)
S_H1 = ParagraphStyle('h1', fontName='DVB', fontSize=13, leading=16, spaceAfter=2)
S_SUB = ParagraphStyle('s', fontName='DV', fontSize=8.5, leading=11, textColor=colors.HexColor('#555555'))
S_H2 = ParagraphStyle('h2', fontName='DVB', fontSize=10.5, leading=13, spaceBefore=7, spaceAfter=3)
S_CAP = ParagraphStyle('c', fontName='DV', fontSize=8.2, leading=10.5, textColor=colors.HexColor('#333333'),
                       spaceBefore=2, spaceAfter=6)
S_CELL = ParagraphStyle('t', fontName='DV', fontSize=7.6, leading=9.2)
S_BUL = ParagraphStyle('bl', parent=S_BODY, leftIndent=10, bulletIndent=0)

D = json.loads((RES / 'fe_results.json').read_text())
FORM = D['formula']
ROWS = D['rows']
NAMES = ['6306', '6210', '6310']


def vn(x, nd=2):
    s = f'{x:,.{nd}f}'.replace(',', ' ').replace('.', ',')
    return s.replace('-', '−')


def pct(fe, f):
    return 100.0 * (f - fe) / fe


def vp(x, nd=1):
    s = f'{x:+.{nd}f}'.replace('.', ',').replace('-', '−')
    return s


def row(name, case, shoulder=0.2):
    for r in ROWS:
        if r['bearing'] == name and r['case'] == case and abs(r['shoulder'] - shoulder) < 1e-9:
            return r
    raise KeyError((name, case, shoulder))


def P(t, st=S_BODY):
    return Paragraph(t, st)


S_HEAD = ParagraphStyle('th', fontName='DVB', fontSize=7.4, leading=8.8, alignment=2)
S_HEAD_L = ParagraphStyle('thl', parent=S_HEAD, alignment=0)


def table(data, widths, header_rows=1, align_from=1, zebra=False):
    data = [list(r) for r in data]
    for i in range(header_rows):
        data[i] = [Paragraph(c, S_HEAD_L if j < align_from else S_HEAD) if isinstance(c, str) and c else c
                   for j, c in enumerate(data[i])]
    t = Table(data, colWidths=widths, repeatRows=header_rows)
    st = [('FONT', (0, 0), (-1, -1), 'DV', 7.6),
          ('FONT', (0, 0), (-1, header_rows - 1), 'DVB', 7.6),
          ('LINEABOVE', (0, 0), (-1, 0), 0.8, colors.black),
          ('LINEBELOW', (0, header_rows - 1), (-1, header_rows - 1), 0.5, colors.black),
          ('LINEBELOW', (0, -1), (-1, -1), 0.8, colors.black),
          ('ALIGN', (align_from, 0), (-1, -1), 'RIGHT'),
          ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
          ('TOPPADDING', (0, 0), (-1, -1), 1.6), ('BOTTOMPADDING', (0, 0), (-1, -1), 1.6),
          ('LEFTPADDING', (0, 0), (-1, -1), 2.5), ('RIGHTPADDING', (0, 0), (-1, -1), 2.5)]
    t.setStyle(TableStyle(st))
    return t


# --------------------------------------------------------------------- derived numbers
QTY = [('p_f', 'p_f', 'p_f, MPa', 2), ('s_bore', 's_bore_mid', 'σθ,bore, MPa', 2),
       ('s_z0', 's_z0', 'σθ,z0, MPa', 2), ('k_cl', 'k_cl', 'k_cl', 3)]
dev = {n: {q: (pct(row(n, 'C')[kr], FORM[n]['ref9'][kf]), pct(row(n, 'C')[kr], FORM[n]['equal_area'][kf]))
           for kf, kr, q, _ in QTY} for n in NAMES}


def rng(q, i):
    v = [dev[n][q][i] for n in NAMES]
    return min(v), max(v)


p9 = rng('p_f, MPa', 0)
pea = rng('p_f, MPa', 1)
sb9 = rng('σθ,bore, MPa', 0)
sbea = rng('σθ,bore, MPa', 1)
sz9 = rng('σθ,z0, MPa', 0)
szea = rng('σθ,z0, MPa', 1)
k9 = rng('k_cl', 0)
kea = rng('k_cl', 1)
ratio_p = [row(n, 'C')['p_f'] / FORM[n]['ref9']['p_f'] for n in NAMES]
ratio_z = {n: row(n, 'C')['s_z0'] / FORM[n]['ref9']['s_z0'] for n in NAMES}
dl0 = delta_max_life(B6310, 5200.0, 18000.0, 1500.0)
dl1 = 4.0 + (dl0 - 4.0) / ratio_z['6310']
caseA = {n: row(n, 'A') for n in NAMES}
errA = max(max(abs(pct(FORM[n]['ref9']['p_f'], caseA[n]['p_f'])),
               abs(pct(FORM[n]['ref9']['s_bore'], caseA[n]['s_bore_mid'])),
               abs(pct(FORM[n]['ref9']['k_cl'], caseA[n]['k_cl']))) for n in NAMES)
conv = D['convergence']
lin = D['linearity_30um']
ext = D['shaft_ext_2B']
c6310 = row('6310', 'C')
shoulder_dz = max(abs(pct(row(n, 'C')['s_z0'], row(n, 'D', s)['s_z0'])) for n in NAMES for s in (0.15, 0.25))
shoulder_dp = max(abs(pct(row(n, 'C')['p_f'], row(n, 'D', s)['p_f'])) for n in NAMES for s in (0.15, 0.25))
shoulder_dk = max(abs(pct(row(n, 'C')['k_cl'], row(n, 'D', s)['k_cl'])) for n in NAMES for s in (0.15, 0.25))
p9_d015 = [pct(row(n, 'D', 0.15)['p_f'], FORM[n]['ref9']['p_f']) for n in NAMES]
p9_d025 = [pct(row(n, 'D', 0.25)['p_f'], FORM[n]['ref9']['p_f']) for n in NAMES]
all_comp = all(r['all_compressive'] for r in ROWS)


def build(path):
    doc = SimpleDocTemplate(str(path), pagesize=A4, leftMargin=17 * mm, rightMargin=17 * mm,
                            topMargin=15 * mm, bottomMargin=15 * mm,
                            title='Kiểm nhánh lắp ghép bằng phần tử hữu hạn', author='')
    W = A4[0] - 34 * mm
    s = []
    s.append(P('Kiểm nhánh lắp ghép bằng phần tử hữu hạn (Hướng 1, Bài báo 2, việc 3)', S_H1))
    s.append(P('13/11/2026 · Ngoài đề cương luận văn: chỉ dùng cho Bài báo 2, <b>chờ thầy duyệt</b>. '
               'Mô hình đối xứng trục đàn hồi tuyến tính, mã nguồn mở (scikit-fem 12.0.2, triangle).', S_SUB))
    s.append(Spacer(1, 4))

    # ------------------------------------------------------------- 1. key results
    s.append(P('1. Kết quả chính', S_H2))
    bl = [
        f'Ca A (vòng trơn, đường kính ngoài = D<sub>eff</sub> theo [9]): FE trùng nghiệm Lamé, sai lệch lớn nhất '
        f'{vn(errA, 4)} % cho p<sub>f</sub>, σθ,bore, k<sub>cl</sub>. Mã FE đạt tiêu chí ≤ 0,5 %.',
        f'Ca C, D<sub>eff</sub> theo [9]: σθ,bore lệch {vp(sb9[0])} … {vp(sb9[1])} %, σθ,z0 lệch '
        f'{vp(sz9[0])} … {vp(sz9[1])} %, k<sub>cl</sub> = d/D<sub>eff</sub> lệch {vp(k9[0])} … {vp(k9[1])} %. '
        f'Ba đại lượng này <b>khớp</b> (≤ 5 %) ở cả ba ổ.',
        f'Ca C, p<sub>f</sub> trung bình: công thức với D<sub>eff</sub> theo [9] thấp hơn FE '
        f'{vn(-p9[1], 1)} … {vn(-p9[0], 1)} %: <b>không khớp</b> (&gt; 10 %). Tỉ số FE/công thức '
        f'{vn(min(ratio_p), 3)} … {vn(max(ratio_p), 3)}.',
        f'D<sub>eff</sub> cùng diện tích cho p<sub>f</sub> lệch {vp(pea[0])} … {vp(pea[1])} % (khớp) nhưng σθ,bore '
        f'{vp(sbea[0])} … {vp(sbea[1])} %, σθ,z0 {vp(szea[0])} … {vp(szea[1])} %, k<sub>cl</sub> '
        f'{vp(kea[0])} … {vp(kea[1])} %. Không có một D<sub>eff</sub> nào khớp cả bốn đại lượng: D<sub>eff</sub> '
        f'theo [9] gần FE hơn ở σθ và k<sub>cl</sub>; D<sub>eff</sub> cùng diện tích gần FE hơn ở p<sub>f</sub>.',
        f'Chiều sai lệch với tuổi thọ: công thức cho σθ,z0 thấp hơn FE {vn(-sz9[1], 1)} … {vn(-sz9[0], 1)} %, tức '
        f'hơi thiên về không an toàn. Ở ca tham chiếu 6310, F<sub>r</sub> = 5 200 N, L<sub>req</sub> = 18 000 h, '
        f'nhân K<sub>δ</sub> với tỉ số FE/công thức ({vn(ratio_z["6310"], 4)}) thì δ<sub>life</sub> đổi từ '
        f'{vn(dl0)} µm thành {vn(dl1)} µm.',
        f'p<sub>f</sub> không đi thẳng vào chuỗi tuổi thọ: chương trình chỉ dùng tích p<sub>f</sub>·g(D<sub>eff</sub>) '
        f'= σθ, còn δ<sub>slip</sub> lấy theo công thức của NTN [4], không từ p<sub>f</sub>. Sai lệch p<sub>f</sub> '
        f'chỉ ảnh hưởng tới các giá trị p<sub>f</sub> được in ra trong bảng.',
    ]
    for b_ in bl:
        s.append(Paragraph(b_, S_BUL, bulletText='•'))

    # ------------------------------------------------------------- 2. setup
    s.append(P('2. Thiết lập', S_H2))
    s.append(P(
        'Định nghĩa D<sub>eff</sub> trong [9] (NASA/TM-2010-216913): tr. 4, pt. (6) – (8a); Hình 6, tr. 5; ví dụ '
        'tính cho ổ cỡ 210, bước (2), tr. 16: D<sub>eff</sub> = 58,94 mm. Theo văn bản, D<sub>eff</sub> được tính '
        '«bằng cách cộng diện tích vai vào vòng không vai rồi chia cho bề rộng». Hình 6 ghi chiều cao vai 0,2·d '
        'theo bán kính. Pt. (8a) cộng (0,2·d − A<sub>sector</sub>/W) một lần vào <i>đường kính</i>, nên phần tăng '
        'bán kính chỉ bằng một nửa của vòng trơn cùng diện tích. Vòng trơn cùng diện tích mặt cắt có '
        'D = D<sub>IR</sub> + 0,4·D<sub>w</sub> − 2A/B (A = 0,1118·D<sub>w</sub>² như chương trình). Gói '
        '<font name="DV">bearing_fit</font> được dựng lại từ bản Drive 27/10 và đạt hai phép kiểm hồi quy '
        '(f<sub>fit</sub> = 0,8591; δ<sub>max,life</sub> = 15,70 µm).'))
    s.append(P(
        'Mô hình FE: đối xứng trục, nửa z ≥ 0, phần tử tam giác bậc hai, E = 205 878 MPa, ν = 0,30. Trục đặc '
        'cùng vật liệu. Lưới trùng nút trên bề mặt lắp. Ràng buộc u<sub>r</sub>(vòng) − u<sub>r</sub>(trục) = '
        'δ<sub>eff</sub>/2 tại mọi bậc tự do hướng kính, áp bằng nhân tử Lagrange; u<sub>z</sub> hai bên độc lập '
        '(không ma sát). Biên: u<sub>r</sub> = 0 tại r = 0, u<sub>z</sub> = 0 tại z = 0. '
        'p<sub>f</sub> = Σλ / [(d/2)(B/2)]. σθ chiếu L2 về trường P2 liên tục rồi lấy tại điểm. Vòng có rãnh: '
        'bán kính cung 0,52·D<sub>w</sub>, đáy rãnh tại D<sub>IR</sub>/2, vai sâu 0,2·D<sub>w</sub> (giả định), '
        'không mép vát. Ca A: vòng trơn D<sub>eff</sub> theo [9], trục dài bằng vòng. Ca B: rãnh thật, trục dài '
        'bằng vòng. Ca C: trục dài thêm B mỗi phía. Ca D: như C, vai 0,15 và 0,25·D<sub>w</sub>. '
        'δ<sub>eff</sub> = 10 µm. D<sub>z</sub> từ <font name="DV">closed_form_constants</font> tại '
        'F<sub>r</sub>/C = 0,08. Sai lệch = (công thức − FE)/FE.'))
    tick = 'Có' if all_comp else 'Không'
    s.append(P(f'Tiếp xúc kín được kiểm sau khi giải: mọi nhân tử đều nén ở mọi ca ({tick}).'))

    # ------------------------------------------------------------- 3. main table
    s.append(P('3. Bảng chính: ca C, δ<sub>eff</sub> = 10 µm', S_H2))
    hdr1 = ['Ổ'] + sum([[q, '', ''] for _, _, q, _ in QTY], [])
    hdr2 = [''] + ['FE', '[9]', 'cùng d.t.'] * 4
    data = [hdr1, hdr2]
    for n in NAMES:
        r = row(n, 'C')
        line = [n]
        for kf, kr, q, nd in QTY:
            line += [vn(r[kr], nd),
                     P(f'{vn(FORM[n]["ref9"][kf], nd)}<br/>({vp(dev[n][q][0])})', S_CELL),
                     P(f'{vn(FORM[n]["equal_area"][kf], nd)}<br/>({vp(dev[n][q][1])})', S_CELL)]
        data.append(line)
    w0 = 11 * mm
    t = table(data, [w0] + [(W - w0) / 12] * 12, header_rows=2)
    t.setStyle(TableStyle([('SPAN', (1, 0), (3, 0)), ('SPAN', (4, 0), (6, 0)), ('SPAN', (7, 0), (9, 0)),
                           ('SPAN', (10, 0), (12, 0)), ('ALIGN', (1, 0), (-1, 1), 'CENTER'),
                           ('LINEBELOW', (1, 0), (-1, 0), 0.4, colors.grey)]))
    s.append(t)
    s.append(P('σθ,bore lấy tại mặt phẳng giữa. k<sub>cl</sub> theo FE = 2·u<sub>r</sub>(đáy rãnh, z = 0)/δ<sub>eff</sub>; '
               'theo công thức = d/D<sub>eff</sub> (mô phỏng 10/11). «[9]»: công thức với D<sub>eff</sub> theo [9]; '
               '«cùng d.t.»: với D<sub>eff</sub> cùng diện tích. Trong ngoặc: sai lệch so với FE, %.', S_CAP))
    data = [['Ổ', 'D_IR', 'D_eff [9]', 'D_eff cùng d.t.', 'D cùng d.t. (mặt cắt FE)', 'D tính ngược từ p_f FE',
             'σθ,bore FE trung bình bề rộng']]
    for n in NAMES:
        r = row(n, 'C')
        data.append([n, vn(FORM[n]['D_IR']), vn(FORM[n]['ref9']['D_eff']), vn(FORM[n]['equal_area']['D_eff']),
                     vn(r['D_eq_area_exact']), vn(r['D_eff_backcalc_p']),
                     f'{vn(r["s_bore_avg"])} ({vp(pct(r["s_bore_avg"], FORM[n]["ref9"]["s_bore"]))} %)'])
    s.append(table([['Ổ', 'D<sub>IR</sub>', 'D<sub>eff</sub> [9]', 'D<sub>eff</sub> cùng d.t.',
                     'D cùng d.t. (mặt cắt FE)', 'D tính ngược từ p<sub>f</sub> FE',
                     'σθ,bore FE trung bình bề rộng']] + data[1:],
                   [11 * mm, 17 * mm, 20 * mm, 24 * mm, 30 * mm, 30 * mm, W - 132 * mm]))
    s.append(P('Đường kính tính bằng mm. «Mặt cắt FE»: diện tích đúng của rãnh bán kính 0,52·D<sub>w</sub>. '
               '«Tính ngược»: D làm pt. (3.9) cho đúng p<sub>f</sub> FE. Cột cuối: σθ trung bình dọc lỗ, MPa, '
               'trong ngoặc là sai lệch của công thức [9].', S_CAP))

    # ------------------------------------------------------------- figures
    img = Image(str(RES / 'fig_field_6310C.png'), width=W * 0.86, height=W * 0.86 * 1061 / 1561)
    s.append(KeepTogether([img, P('Hình 1. Trường σθ trong vòng trong 6310, ca C, δ<sub>eff</sub> = 10 µm, nửa '
                                  'trên mặt phẳng giữa. Ba điểm lấy kết quả: 1 lỗ; 2 độ sâu z<sub>0</sub>; '
                                  '3 đáy rãnh.', S_CAP)]))
    img = Image(str(RES / 'fig_radial_6310C.png'), width=W * 0.78, height=W * 0.78 * 0.6)
    s.append(KeepTogether([img, P('Hình 2. σθ dọc bán kính tại z = 0, 6310, ca C, so với nghiệm Lamé với hai '
                                  'giá trị D<sub>eff</sub>.', S_CAP)]))

    # ------------------------------------------------------------- 4. supplementary tables
    s.append(P('4. Bảng phụ', S_H2))
    data = [['Ổ', 'p<sub>f</sub> FE', 'p<sub>f</sub> Lamé', 'σθ,bore FE', 'σθ,bore Lamé',
             'k<sub>cl</sub> FE (mặt ngoài)', 'd/D<sub>eff</sub>', 'sai lệch lớn nhất']]
    for n in NAMES:
        r, f = caseA[n], FORM[n]['ref9']
        e = max(abs(pct(f['p_f'], r['p_f'])), abs(pct(f['s_bore'], r['s_bore_mid'])),
                abs(pct(f['k_cl'], r['k_cl'])))
        data.append([n, vn(r['p_f'], 4), vn(f['p_f'], 4), vn(r['s_bore_mid'], 4), vn(f['s_bore'], 4),
                     vn(r['k_cl'], 5), vn(f['k_cl'], 5), f'{vn(e, 4)} %'])
    s.append(KeepTogether([P('Bảng 3. Ca A, kiểm mã.', S_CAP),
                           table(data, [11 * mm] + [(W - 11 * mm) / 7] * 7)]))
    s.append(Spacer(1, 5))
    data = [['Mức', 'cỡ phần tử dưới rãnh, mm', 'cỡ phần tử chung, mm', 'số bậc tự do', 'p<sub>f</sub>',
             'σθ,bore', 'σθ,z0', 'k<sub>cl</sub>', 'Δσθ,z0 so với mức trước']]
    prev = None
    for i, r in enumerate(conv):
        dz = '' if prev is None else f'{vn(abs(pct(prev, r["s_z0"])), 4)} %'
        data.append([str(i + 1), vn(r['h_fine'], 4), vn(r['h_glob'], 3), f'{r["n_dof"]:,}'.replace(',', ' '),
                     vn(r['p_f'], 4), vn(r['s_bore_mid'], 4), vn(r['s_z0'], 4), vn(r['k_cl'], 5), dz])
        prev = r['s_z0']
    z0 = conv[0]['z0']
    s.append(KeepTogether([P(f'Bảng 4. Hội tụ lưới, 6310 ca C (z<sub>0</sub> = {vn(z0, 4)} mm). Mức 3 dùng cho '
                             f'mọi ca; cỡ phần tử dưới rãnh = z<sub>0</sub>/8.', S_CAP),
                           table(data, [10 * mm, 22 * mm, 20 * mm, 20 * mm, 18 * mm, 18 * mm, 18 * mm, 16 * mm,
                                        W - 142 * mm])]))
    s.append(P(
        f'Tuyến tính (6310 ca C, δ<sub>eff</sub> = 30 µm so với 10 µm): tỉ số p<sub>f</sub> = '
        f'{vn(lin["p_f"] / c6310["p_f"], 5)}, σθ,bore = {vn(lin["s_bore_mid"] / c6310["s_bore_mid"], 5)}, '
        f'σθ,z0 = {vn(lin["s_z0"] / c6310["s_z0"], 5)}; k<sub>cl</sub> = {vn(lin["k_cl"], 5)}. '
        f'Trục dài thêm 2B mỗi phía thay vì B: p<sub>f</sub> đổi {vp(pct(c6310["p_f"], ext["p_f"]), 2)} %, '
        f'σθ,z0 đổi {vp(pct(c6310["s_z0"], ext["s_z0"]), 2)} %.', S_CAP))

    data = [['Ổ', 'Ca', 'vai', 'p<sub>f</sub> FE', 'sai lệch [9]', 'σθ,bore FE', 'sai lệch [9]', 'σθ,z0 FE',
             'sai lệch [9]', 'k<sub>cl</sub> FE', 'sai lệch [9]']]
    for n in NAMES:
        f = FORM[n]['ref9']
        for case, sh in (('B', 0.2), ('C', 0.2), ('D', 0.15), ('D', 0.25)):
            r = row(n, case, sh)
            data.append([n, case, vn(sh, 2), vn(r['p_f'], 3), f'{vp(pct(r["p_f"], f["p_f"]))} %',
                         vn(r['s_bore_mid'], 2), f'{vp(pct(r["s_bore_mid"], f["s_bore"]))} %',
                         vn(r['s_z0'], 2), f'{vp(pct(r["s_z0"], f["s_z0"]))} %',
                         vn(r['k_cl'], 4), f'{vp(pct(r["k_cl"], f["k_cl"]))} %'])
    t = table(data, [10 * mm, 8 * mm, 10 * mm] + [(W - 28 * mm) / 8] * 8, align_from=2)
    t.setStyle(TableStyle([('LINEBELOW', (0, 4), (-1, 4), 0.3, colors.grey),
                           ('LINEBELOW', (0, 8), (-1, 8), 0.3, colors.grey)]))
    s.append(KeepTogether([P('Bảng 5. Ca B, C, D: giá trị FE và sai lệch của công thức với D<sub>eff</sub> theo '
                             '[9] (công thức không đổi theo chiều sâu vai). Vai tính theo D<sub>w</sub>.', S_CAP), t]))
    s.append(P(
        f'Thăm dò vai (ca D, ±0,05·D<sub>w</sub>): σθ,z0 đổi tối đa {vn(shoulder_dz, 2)} %, k<sub>cl</sub> '
        f'{vn(shoulder_dk, 2)} %, p<sub>f</sub> {vn(shoulder_dp, 1)} %. Sai lệch p<sub>f</sub> của [9] là '
        f'{vp(min(p9_d015))} … {vp(max(p9_d015))} % với vai 0,15·D<sub>w</sub> và {vp(min(p9_d025))} … '
        f'{vp(max(p9_d025))} % với vai 0,25·D<sub>w</sub>. Kết luận về σθ và k<sub>cl</sub> không phụ thuộc giả '
        f'định vai; mức sai lệch p<sub>f</sub> thì phụ thuộc.', S_CAP))
    img = Image(str(RES / 'fig_pressure_6310.png'), width=W * 0.72, height=W * 0.72 * 0.49)
    s.append(KeepTogether([img, P('Hình 3. Áp suất tiếp xúc nút dọc bề mặt lắp, 6310, chuẩn hoá theo trung bình. '
                                  'Ca C có tập trung áp suất ở mép vòng; đỉnh phụ thuộc lưới nên chỉ dùng giá trị '
                                  'trung bình.', S_CAP)]))

    # ------------------------------------------------------------- 5. limits
    s.append(P('5. Giới hạn', S_H2))
    for b_ in [
        'Biên dạng vai và chiều sâu vai 0,2·D<sub>w</sub> là giả định, không có bản vẽ của hãng; bỏ qua mép vát. '
        f'Ca D cho thấy p<sub>f</sub> nhạy với giả định này (tới {vn(shoulder_dp, 1)} %), σθ và k<sub>cl</sub> '
        f'thì không (≤ {vn(max(shoulder_dz, shoulder_dk), 2)} %).',
        'Bề mặt lắp không ma sát, tiếp xúc kín toàn bề rộng (đã kiểm: mọi nhân tử nén), trục đặc cùng vật liệu, '
        'đàn hồi tuyến tính, biến dạng nhỏ, không nhiệt, không tải bi.',
        'Áp suất ở mép bề mặt lắp (ca C) là điểm kỳ dị của mô hình đàn hồi: đỉnh tăng khi làm mịn lưới. '
        'p<sub>f</sub> lấy trung bình từ tổng nhân tử, đã hội tụ (Bảng 4).',
        'FE chỉ kiểm ứng suất và biến dạng do lắp ghép. Không kiểm chứng tuổi thọ, không thay được thực nghiệm.',
    ]:
        s.append(Paragraph(b_, S_BUL, bulletText='•'))

    # ------------------------------------------------------------- 6. open points
    s.append(P('6. Việc chờ tác giả và thầy', S_H2))
    for b_ in [
        'Thầy có đồng ý đưa mục kiểm FE vào Bài báo 2 không.',
        f'p<sub>f</sub> theo D<sub>eff</sub> của [9] lệch FE trên 10 % ({vn(-p9[1], 1)} … {vn(-p9[0], 1)} %). '
        'Bài báo 1 và luận văn (ví dụ Bảng 4.3) có in giá trị p<sub>f</sub> không? Nếu có, tác giả quyết định cách '
        'xử lý. Báo cáo này không sửa gì trong chương trình hay bản thảo.',
        'Tác giả có bản vẽ hoặc số đo vai, mép vát của vòng trong 6306, 6210, 6310 không? Nếu có, chạy lại ca C '
        'với hình học đó.',
    ]:
        s.append(Paragraph(b_, S_BUL, bulletText='•'))
    s.append(P('Tài liệu: [4] NTN, catalog ổ lăn; [9] NASA/TM-2010-216913, Interference-fit life factors for ball '
               'bearings, 2010.', S_CAP))

    def footer(canv, doc_):
        canv.saveState()
        canv.setFont('DV', 7.5)
        canv.setFillColor(colors.HexColor('#666666'))
        canv.drawRightString(A4[0] - 17 * mm, 9 * mm, f'{doc_.page}')
        canv.drawString(17 * mm, 9 * mm, 'Kiểm FE nhánh lắp ghép · Bài báo 2 · chờ thầy duyệt')
        canv.restoreState()

    doc.build(s, onFirstPage=footer, onLaterPages=footer)


if __name__ == '__main__':
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else RES / 'BaoCao_KiemFE_Huong1_20261113.pdf'
    build(out)
    print(out)
